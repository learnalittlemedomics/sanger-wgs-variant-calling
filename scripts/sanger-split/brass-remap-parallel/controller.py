#!/usr/bin/env python3
"""Ordered, restartable orchestration of the unchanged BRASS remap analysis."""
import concurrent.futures, fcntl, hashlib, json, os, pathlib, signal, subprocess, sys, threading, time

HERE = pathlib.Path(__file__).resolve().parent
WORKER = HERE / 'worker.pl'
ARRAYS = ['low_end_remap_score','high_end_remap_score','low_end_footprint','high_end_footprint']
CHILDREN = {}; MUTEX = threading.Lock(); STOP = threading.Event()

def digest(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024), b''): h.update(block)
    return h.hexdigest()

def atomic(path, obj):
    path = pathlib.Path(path); temp = path.with_name(path.name+'.pending')
    with temp.open('w') as f:
        json.dump(obj,f,sort_keys=True); f.flush(); os.fsync(f.fileno())
    os.replace(str(temp),str(path))

def stat_identity(path):
    p=pathlib.Path(path).resolve(); s=p.stat()
    return [str(p),s.st_size,s.st_mtime_ns]

def options(argv):
    values={}; pos=[]; i=0
    allowed=set('blat_mapping_threshold groups_file slop_for_genomic_region slop_for_getting_reads small_indel_threshold remapping_score_threshold microbial_fraction_threshold footprint_size_max scores_output_file virus_db bacterial_db_stub tmpdir score_alg search_cores min_support'.split())
    while i<len(argv):
        a=argv[i]
        if a.startswith('-'):
            k=a.lstrip('-')
            if k=='noclean': values[k]=True; i+=1; continue
            if k not in allowed or i+1>=len(argv): raise ValueError('Unsupported option '+a)
            values[k]=argv[i+1];i+=2
        else:pos.append(a);i+=1
    if len(pos)!=4: raise ValueError('Expected BEDPE, BAM, reference and output paths')
    for k in ['groups_file','scores_output_file','tmpdir','virus_db','bacterial_db_stub']:
        if not values.get(k):raise ValueError('Missing '+k)
    if values.get('score_alg','ssearch36')!='ssearch36':raise ValueError('Validated backend is ssearch36')
    return values,pos

def invoke(cmd, directory):
    with (directory/'stdout.log').open('w') as out,(directory/'stderr.log').open('w') as err:
        with MUTEX:
            if STOP.is_set():raise RuntimeError('Worker launch canceled')
            p=subprocess.Popen(cmd,stdout=out,stderr=err,start_new_session=True);CHILDREN[p.pid]=p
        try:
            rc=p.wait()
            if rc:raise RuntimeError('Worker failed rc={} log={}'.format(rc,directory/'stderr.log'))
        finally:
            with MUTEX:CHILDREN.pop(p.pid,None)

def terminate_children():
    STOP.set()
    with MUTEX:
        for pid,p in list(CHILDREN.items()):
            if p.poll() is None:
                try:os.killpg(pid,signal.SIGTERM)
                except ProcessLookupError:pass

def interrupted(signum, frame):
    terminate_children();raise RuntimeError('Interrupted by signal '+str(signum))

def worker_command(opts,pos,d,extra,threads):
    updated=dict(opts);updated.update(groups_file=str(d/'groups.tsv'),tmpdir=str(d/'tmp'),search_cores=str(threads),noclean=True,scores_output_file=str(d/'scores.pending'))
    cmd=['perl',str(WORKER)]
    for k,v in sorted(updated.items()):
        cmd.append('-'+k)
        if v is not True:cmd.append(str(v))
    return cmd+extra+[str(d/'input.r4'),pos[1],pos[2],str(d/'result.pending')]

def verify_chunk(d, expected):
    receipt=json.loads((d/'complete.json').read_text())
    if receipt['input_sha256']!=expected['input_sha256'] or receipt['groups_sha256']!=expected['groups_sha256']:raise RuntimeError('Chunk identity mismatch '+str(d))
    for name,wanted in receipt['files'].items():
        p=(d/name).resolve()
        if d.resolve() not in p.parents or digest(p)!=wanted:raise RuntimeError('Corrupt checkpoint '+str(p))
    state=json.loads((d/'state.json').read_text())
    if any(len(state[k])!=expected['records'] for k in ARRAYS):raise RuntimeError('Chunk count mismatch')
    if len(state['all_reads_files'])!=1:raise RuntimeError('Unexpected reads list')
    p=pathlib.Path(state['all_reads_files'][0]).resolve()
    if d.resolve() not in p.parents or str(p.relative_to(d.resolve())) not in receipt['files']:raise RuntimeError('Unverified reads file')
    return state

def main(argv):
    opts,pos=options(argv)
    cores=int(opts.get('search_cores','1'));workers=min(15,max(1,cores//2));threads=min(2,cores)
    if cores<1:raise ValueError('Invalid CPU allocation')
    # The testing override only changes partition boundaries, never scoring.
    chunk_size=int(os.environ.get('BRASS_REMAP_CHUNK_SIZE','256'))
    if chunk_size<1:raise ValueError('Invalid chunk size')
    root=pathlib.Path(opts['tmpdir']).absolute()/'parallel-v1';root.mkdir(parents=True,exist_ok=True)
    lock=(root/'controller.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    inputs={'bedpe_sha256':digest(pos[0]),'groups_sha256':digest(opts['groups_file']), 'bam':stat_identity(pos[1]),'reference':stat_identity(pos[2]),'argv':argv,'worker_sha256':digest(WORKER),'controller_sha256':digest(__file__),'chunk_size':chunk_size}
    inputs['databases']=[stat_identity(opts['virus_db'])]+[stat_identity(opts['bacterial_db_stub']+'.'+str(i)+'.fa.2bit') for i in range(1,5)]
    identity=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    if (root/'manifest.json').exists():
        manifest=json.loads((root/'manifest.json').read_text())
        if manifest['identity']!=identity:raise RuntimeError('Checkpoint belongs to different inputs/code/options')
    else:
        chunks=[];id_chunks={};current=[]
        def save_rows(rows):
            number=len(chunks);d=root/('chunk-{:05d}'.format(number));d.mkdir(exist_ok=True)
            (d/'tmp').mkdir(exist_ok=True);(d/'input.r4').write_text(''.join(rows))
            for line in rows:id_chunks.setdefault(line.rstrip('\n').split('\t')[6],set()).add(number)
            chunks.append({'name':d.name,'records':len(rows),'input_sha256':digest(d/'input.r4')})
        with open(pos[0]) as f:
            for line in f:
                if line.startswith('#'):continue
                if not line.endswith('\n'):line+='\n'
                current.append(line)
                if len(current)==chunk_size:save_rows(current);current=[]
        if current or not chunks:save_rows(current)
        # Append each groups row to every chunk needing it, preserving file order.
        for c in chunks:(root/c['name']/'groups.tsv').write_text('')
        handles={}
        try:
            with open(opts['groups_file']) as f:
                for line in f:
                    if line.startswith('#'):continue
                    group_id=line.rstrip('\n').split('\t')[6]
                    for number in sorted(id_chunks.get(group_id,())):
                        if number not in handles:
                            if len(handles)>=64:
                                _,handle=handles.popitem();handle.close()
                            handles[number]=(root/chunks[number]['name']/'groups.tsv').open('a')
                        handles[number].write(line)
        finally:
            for handle in handles.values():handle.close()
        for c in chunks:c['groups_sha256']=digest(root/c['name']/'groups.tsv')
        manifest={'identity':identity,'inputs':inputs,'chunks':chunks,'total_records':sum(c['records'] for c in chunks)}
        atomic(root/'manifest.json',manifest)
    if (root/'complete.json').exists():
        done=json.loads((root/'complete.json').read_text())
        if done['identity']!=identity or any(not pathlib.Path(p).is_file() or digest(p)!=h for p,h in done['outputs'].items()):raise RuntimeError('Final checkpoint/output mismatch')
        print('REMAP_PARALLEL_COMPLETE_CACHED records={}'.format(done['records']),flush=True);return
    chunks=manifest['chunks'];completed=set();completed_records=0
    for c in chunks:
        d=root/c['name']
        if digest(d/'input.r4')!=c['input_sha256'] or digest(d/'groups.tsv')!=c['groups_sha256']:raise RuntimeError('Changed chunk inputs')
        if (d/'complete.json').exists():verify_chunk(d,c);completed.add(c['name']);completed_records+=c['records']
    print('REMAP_PARALLEL workers={} threads={} records={} cached_chunks={}'.format(workers,threads,manifest['total_records'],len(completed)),flush=True)
    def progress(stage):
        atomic(root/'progress.json',{'stage':stage,'completed_chunks':len(completed),'total_chunks':len(chunks),'completed_records':completed_records,'total_records':manifest['total_records'],'workers':workers,'threads_per_worker':threads,'timestamp':time.time()})
    progress('remapping')
    def run(c):
        d=root/c['name']
        extra=['-remap_state',str(d/'state.pending.json')]
        invoke(worker_command(opts,pos,d,extra,threads),d)
        os.replace(str(d/'state.pending.json'),str(d/'state.json'))
        state=json.loads((d/'state.json').read_text())
        if any(len(state[k])!=c['records'] for k in ARRAYS):raise RuntimeError('Worker incomplete arrays')
        files={'state.json':digest(d/'state.json')}
        for p in state['all_reads_files']:
            p=pathlib.Path(p).resolve()
            if d.resolve() not in p.parents:raise RuntimeError('Worker reads outside chunk')
            files[str(p.relative_to(d.resolve()))]=digest(p)
        atomic(d/'complete.json',{'input_sha256':c['input_sha256'],'groups_sha256':c['groups_sha256'],'files':files,'records':c['records']})
        return c
    pool=concurrent.futures.ThreadPoolExecutor(max_workers=workers)
    try:
        futures=[pool.submit(run,c) for c in chunks if c['name'] not in completed]
        for future in concurrent.futures.as_completed(futures):
            c=future.result();completed.add(c['name']);completed_records+=c['records'];progress('remapping')
    except BaseException:
        terminate_children();progress('failed');raise
    finally:pool.shutdown(wait=True)
    merged={k:[] for k in ARRAYS+['all_reads_files']};merged['rg_id_of_read']={}
    for c in chunks:
        state=verify_chunk(root/c['name'],c)
        for k in ARRAYS+['all_reads_files']:merged[k].extend(state[k])
        merged['rg_id_of_read'].update(state['rg_id_of_read'])
    atomic(root/'merged-state.json',merged);del merged
    progress('microbial_filter')
    final=root/'final';final.mkdir(exist_ok=True);(final/'tmp').mkdir(exist_ok=True)
    # Full original input and group file, exact original order and parameters.
    updated=dict(opts);updated.update(tmpdir=str(final/'tmp'),scores_output_file=str(final/'scores.pending'))
    cmd=['perl',str(WORKER)]
    for k,v in sorted(updated.items()):
        cmd.append('-'+k)
        if v is not True:cmd.append(str(v))
    cmd+=['-load_remap_state',str(root/'merged-state.json')]+pos[:3]+[str(final/'result.pending')]
    invoke(cmd,final)
    final_files={pos[3]:digest(final/'result.pending'),opts['scores_output_file']:digest(final/'scores.pending')}
    os.replace(str(final/'scores.pending'),opts['scores_output_file'])
    os.replace(str(final/'result.pending'),pos[3])
    atomic(root/'complete.json',{'identity':identity,'outputs':final_files,'records':manifest['total_records']});progress('complete')
    print('REMAP_PARALLEL_COMPLETE records={}'.format(manifest['total_records']),flush=True)

if __name__=='__main__':
    signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
    try:main(sys.argv[1:])
    except BaseException:
        terminate_children();raise
