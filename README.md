# Sanger WGS Variant Calling — ChatGPT 辅助优化版

本仓库是 [icgc-argo-workflows/sanger-wgs-variant-calling](https://github.com/icgc-argo-workflows/sanger-wgs-variant-calling) 的个人 fork，基于上游提交 [817445b](https://github.com/icgc-argo-workflows/sanger-wgs-variant-calling/commit/817445baef7f6579f7240a082a8779cf2fec2b29)（流程版本 `2.1.0-9.7.1`）进行修改。

原项目提供基于 Nextflow 的 Sanger WGS 体细胞变异检测流程，包含 CaVEMan、Pindel、ASCAT 和 BRASS。底层容器化分析工具源自 [cancerit/dockstore-cgpwgs](https://github.com/cancerit/dockstore-cgpwgs)。

本 fork 在原项目基础上，使用 **ChatGPT 辅助分析、编写和整理优化代码**，主要针对流程拆分、并行调度、资源配置、断点续跑和文件处理问题。仓库维护者已在自己的运行环境中执行过修改后的流程。这里的优化针对工作流组织和执行方式，不代表重新实现或改进了底层变异检测算法。

## 具体优化了什么

| 问题或改进目标 | 本 fork 的改动 |
| --- | --- |
| 原 Sanger 调用模块以较大的整体任务执行，阶段之间难以分别调度和复用 | 将准备、参考缓存、ASCAT、Pindel、CaVEMan、BRASS、QC 和结果打包拆分为独立的 Nextflow 任务，使各阶段可分别配置资源，并在满足 Nextflow 缓存条件时通过 `-resume` 复用结果。 |
| CaVEMan、Pindel 和 BRASS 的重计算阶段需要更细粒度的并行控制 | 增加 scatter/gather 执行：CaVEMan 拆分计算阶段并汇总，Pindel 按染色体处理，BRASS 按阶段拆分和汇总；使用分组与受限 worker 并发控制任务规模。 |
| 各阶段的 CPU、内存需求不同，统一分配不够灵活 | 增加按阶段配置的资源表和 scatter 参数；部分失败重试可调整资源或 worker 数。仓库中的数值是可调整的示例，不是适用于所有数据的性能保证。 |
| 部分准备任务等待了不必要的上游步骤 | 调整依赖关系，让满足输入要求的准备和缓存步骤更早启动，同时保留调用器之间必需的依赖。 |
| 修改无关脚本可能影响任务缓存复用 | 按阶段仅传入实际执行所需的脚本，减少无关脚本变更对任务缓存的影响。 |
| 只想运行部分调用器时，原有后处理链不够灵活 | 增加 `params.sanger_steps` 开关、依赖检查及对应后处理分支。CaVEMan 要求启用 ASCAT 和 Pindel；BRASS 和 tumour QC 要求启用 ASCAT。 |
| 集群配额或提交频率限制可能导致任务提交失败 | 增加提交速率、队列规模及针对明确配额拒绝信息的重试配置；将提交拒绝重试与任务失败重试区分开。 |
| 部分模块的 `publishDir.enabled` 使用字符串参数 | 将开关显式转换为布尔值，避免发布开关类型不明确。 |
| 清理任务同时接收多个同名文件时存在暂存冲突 | 为清理输入分配独立的暂存子目录。 |
| payload 的运行标识使用可读运行名称 | 改用 `workflow.sessionId`，统一使用会话 UUID 作为该运行标识。 |
| 运行环境配置与代码耦合，不便共享 | 移除实际服务器队列、账号配额说明、本地容器路径和样本标识；将共享脚本与配置路径按仓库布局整理。 |

另保留已有的 mock 测试与 Nextflow stub 检查，并提供 `scripts/caveman_seff_report.py` 用于整理作业资源使用情况。上述改动记录在本 fork 的 [PR #1](https://github.com/learnalittlemedomics/sanger-wgs-variant-calling/pull/1) 中。

## 配置与使用说明

- Nextflow 版本要求为 **24.10.0 或更新版本**。
- 流程入口：[`sanger-wgs-variant-calling/main.nf`](sanger-wgs-variant-calling/main.nf)。
- 分阶段资源与调用器开关：[`conf/sanger-split.config`](conf/sanger-split.config)。
- scatter 与提交重试参数：[`conf/sanger-scatter.config`](conf/sanger-scatter.config)。
- 更多说明：[`docs/split-workflow.md`](docs/split-workflow.md)。

请在自己的私有配置中设置输入、参考数据、输出目录、容器、执行器和队列等环境参数。若使用本地 Sanger 容器文件，可通过 `params.sanger_split_sif` 指定。请勿将凭据、真实服务器配置或样本信息提交到公开仓库。

## 验证范围

维护者已在自己的环境中运行过修改后的流程。迁移时已有的 CaVEMan、Pindel 和 BRASS mock 测试通过，Python 语法检查通过，并对待发布改动进行了已知私有标识和凭据格式扫描。

mock 测试使用模拟程序和小型假数据，检查任务组织、文件传递、重试及结果汇总等逻辑，不等同于真实样本的生物学结果验证。公开版本经过配置脱敏和相对路径调整，**这些发布整理后的配置与路径未重新进行端到端运行**。目前没有提供统一数据集上的加速倍数、成本降低比例或变异结果一致性基准。

## 致谢

感谢 [ICGC ARGO 原始工作流](https://github.com/icgc-argo-workflows/sanger-wgs-variant-calling) 和 [Cancer IT / Sanger 底层工具](https://github.com/cancerit/dockstore-cgpwgs) 的开发者。本 fork 是使用 ChatGPT 辅助优化的个人版本，不代表上游项目的官方发布。
