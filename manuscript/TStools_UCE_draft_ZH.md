# TStools：一种用于靶向捕获和基因组浅层测序数据的轻量化、证据约束 UCE 恢复工具

> 英文题目建议：**TStools: lightweight and evidence-aware recovery of ultraconserved elements from target-capture and genome-skimming reads**
>
> 草稿范围：仅包括 TStools 核心兼容流程与 UCE 专用流程。软件描述以 GitHub Release v1.6.1 和 PR #12 的合并实现为准；所有未经正式实验确认的内容均以 `【待补：……】` 标记。

**作者：**【待补：作者姓名与顺序】  
**单位：**【待补：作者单位】  
**通讯作者：**【待补：姓名与邮箱】

## 摘要

超保守元件（ultraconserved elements，UCEs）以保守核心提供跨类群定位锚点，并以相邻侧翼提供系统发育信息位点，已广泛用于非模式生物的靶向捕获和基因组浅层测序研究（Faircloth et al., 2012；Smith et al., 2014）。现有流程通常先进行全样本组装，或将 reads 分配至各 locus 后调用通用组装器；这些方案有效但可能带来较高的内存、临时存储和软件依赖成本，而且恢复 locus 的数量或长度并不直接等同于序列正确性和下游可用性（Bossert et al., 2024；Ortiz et al., 2026）。

本文提出 TStools，一种以 Rust 实现的轻量化、参考引导短 reads 恢复工具。其核心流程保留 GeneMiner/GeneMiner2 的参考引导过滤与加权 de Bruijn 图思想（Xie et al., 2024；Yu et al., 2026）。UCE 专用默认 `fast` 模式将广泛招募、完整双端 fragment 保留、方向与最长精确匹配评估及逐 locus 证据选择合并为一次 FASTQ 扫描；v1.6.1 新增的可选 `auto` 模式仅对快速阶段未选中 fragments 的 loci 再扫描一次，并以完整 probe 面板歧义检查、read 局部比对和 contig—probe 证据约束新增结果。TStools 对低深度 locus 保留全部合格证据，对饱和 locus 按参考位置和末端延伸潜力进行有界选择；随后使用参考锚定的加权 de Bruijn 图、双端分支支持和有界前瞻恢复 UCE 核心及 reads 支持的侧翼。可选 rescue 最多执行两轮，只接受具有独立 fragment 和边界跨越证据的延伸，不以参考序列填补无 reads 支持的缺口。

在核心 MainFilter 基准中，Rust 实现与历史 GeneMiner2 基线产生相同的逐 locus reads 输出，同时将峰值内存降低约 50%，运行时间缩短约 5 倍。在 v1.6.1 的开发期回归中，`auto` 相对 `fast` 未丢失已接受 loci，也未改变共有序列；在蜜蜂和河豚的同个体基因组外部参照检查中分别新增 46 和 23 个 loci，其中可外部评分的新增 loci 有 44/45 和 21/21 获得目标 locus 支持。该检查属于强外部参照而非已知真值。【待补：UCE 已知真值模拟的 identity、完整度、indel 和 chimera 结果。】【待补：正式真实 UCE 数据的 locus occupancy、信息位点和 gene-tree congruence。】【待补：Windows 11 与 macOS 普通 16-GB 笔记本的完成率、峰值内存、时间和临时磁盘占用。】这些结果将用于检验 TStools 是否能在不依赖 Python、SPAdes 或服务器级硬件的条件下完成常规 UCE 数据恢复。

TStools 将 UCE 组装表述为一个 target-specific、evidence-bounded 的恢复问题，其目标不是单独最大化 locus 数量，而是在序列正确性、下游可用性和计算成本之间取得更优平衡。

**关键词：** ultraconserved elements；target capture；genome skimming；de Bruijn graph；read recruitment；lightweight bioinformatics

## 1 引言

UCE 的高度保守核心可用于跨较深进化尺度识别同源位点，而共同捕获的侧翼通常包含更多变异，可支持较浅层级的系统发育推断（Faircloth et al., 2012；Smith et al., 2014）。因此，UCE 数据分析的实际目标并非只找回 bait 或保守核心，而是可靠地恢复核心及其有 reads 支持的侧翼。

UCE 和靶向捕获数据已有两类主要处理思路。PHYLUCE 和 UCEasy 等流程通常先用 SPAdes 等通用组装器构建 contigs，再从 contigs 中识别目标 locus（Bankevich et al., 2012；Faircloth, 2016；Ribeiro et al., 2021）。HybPiper 和 SECAPR 则先按目标分配 reads，再进行逐 locus 组装和序列提取（Johnson et al., 2016；Andermann et al., 2018）。Captus 采用全样本 MEGAHIT 组装后提取目标的路线，可统一处理多类测序数据（Li et al., 2015；Ortiz et al., 2026）；组装后的目标识别与提取策略本身也会影响最终恢复结果（Knyshov et al., 2021）。这些路线各有优势，但其 locus 数量、长度、缺失率、信息位点和 gene-tree congruence 之间并不存在单调关系。直接比较 UCE 组装器的研究显示，SPAdes 可恢复较多 loci，而 HybPiper 的结果可能更长、缺失更少并具有更高的 gene-tree congruence；后者是下游一致性指标，不应被表述为已知真值下的“组装准确率”（Bossert et al., 2024）。针对靶向捕获流程的比较也表明，测序深度和工作流结构会共同影响恢复结果与计算成本（Raza et al., 2023）。

Easy353 和 GeneMiner 为低深度测序数据的目标基因快速恢复提供了轻量化基础（Zhang et al., 2022；Xie et al., 2024），GeneMiner2 进一步引入两级哈希过滤、基于方向和结构异常的精细 reads 选择以及自适应 k-mer 组装（Yu et al., 2026）。然而，UCE 与常规编码基因不同：短保守核心可能招募大量重复或多重匹配 reads，而有效的系统发育信号往往位于侧翼；简单限制 reads 数量可能丢失位置多样性和末端延伸证据。另一方面，为每个样本或 locus 启动通用组装器会增加内存、临时文件和依赖管理成本，使普通个人电脑上的批量分析更困难。

为此，我们开发了 TStools。其设计目标是以最少必要机制完成可靠恢复：第一，提供与 GeneMiner2 语义兼容且资源占用更低的核心过滤与组装路径；第二，默认以一次 FASTQ 扫描完成 UCE 招募、证据评估和逐 locus 自适应选择，并仅在显式启用 `auto` 时对未恢复 loci 执行一次受约束的敏感扫描；第三，只允许 reads 明确支持的组装与 rescue 延伸。我们检验三个假设：（i）核心 Rust 路径可在保持输出兼容的同时降低时间和内存；（ii）UCE 专用证据选择、保守敏感招募与组装可提高已知真值下的正确性—完整度平衡；（iii）完整流程可在普通 Windows 和 macOS 笔记本上运行，而不需要 Python、SPAdes 或超大内存服务器。

## 2 材料与方法

### 2.1 软件范围与实现

本文以 TStools v1.6.1（Git commit `3d198e3`）作为软件快照；代码发布于 https://github.com/GUIBA-EX/TStools，版本说明与源代码归档见 https://github.com/GUIBA-EX/TStools/releases/tag/v1.6.1，【待补：永久归档 DOI】。生产代码以 Rust 实现，最低支持 Rust 1.87；【待补：正式实验实际使用的 rustc 完整版本】。本文仅评估两条路径：常规目标的核心兼容路径 `MainFilter → refilter → original-rust`，以及 UCE 专用路径 `UCEFilter fast [→ optional auto fallback] → uce-rust → optional rescue`。除系统压缩库【待补：正式环境中的 zlib 或 zlib-ng 版本】外，生产路径不需要 Python 运行时；UCE 路径不调用 SPAdes 或其他外部通用组装器。

所有分析均记录软件版本、完整命令、输入文件信息、参考与样本表的 SHA-256、关键参数和工作流状态。【待补：正式发布平台，包括 Windows x86-64、macOS Intel/Apple Silicon 和 Linux x86-64；若未完成则删除相应平台。】

### 2.2 核心兼容路径

MainFilter 使用 canonical 双链 2-bit rolling k-mers 招募 reads；任一 mate 命中时保留完整 paired fragment。refilter 根据精确匹配、方向一致性和深度限制进一步选择 reads，随后 `original-rust` 以参考定位的 seeds 构建加权 de Bruijn 图并进行双向延伸。de Bruijn 图组装的基础见 Idury and Waterman（1995）；参考引导的加权节点框架来自 GeneMiner，并以 GeneMiner2 作为算法与结果兼容基线（Xie et al., 2024；Yu et al., 2026）。TStools 的 Rust 实现重写了 k-mer 编码、缓存、FASTX 解析和有界 I/O，但不把这些工程改写表述为 GeneMiner/GeneMiner2 算法的重新发明。

核心兼容性以相同 paired FASTQ、参考、参数和输出模式检验。比较项目包括逐 locus 非空文件数、reads 数、逐文件字节一致性、contig 序列一致性、wall time 和峰值 resident set size（RSS）。【待补：完整 `original-rust` 与历史 GeneMiner2 组装结果对照；当前已完成的实测仅覆盖 MainFilter。】

### 2.3 UCE 专用流程

#### 2.3.1 默认单次扫描的招募与证据评估

默认 `--uce-recruit-mode fast` 使用 *k*=31、step=4 的 UCEFilter 扫描，并以独立的 verification *k*=19 评估精确证据。UCEFilter 预先从所有参考 locus 构建 canonical rolling-k-mer 索引。查询首先经过 cache-local blocked Bloom filter；所有阳性均由精确哈希表复核，因此 Bloom filter 只减少无效查询，不改变招募结果（Bloom, 1970；Putze et al., 2009）。一个 fragment 命中多个 loci 时保留候选集合，但 fragment 本身仅存储一次，R1 与 R2 始终以不可分割单元被选择或丢弃。v1.6.1 将 coarse recruitment *k* 与 verification *k* 解耦，同时保持 `fast` 的默认参数和既有行为不变。

每个 locus 另建包含正向和反向互补序列的 FM-index，以 Burrows-Wheeler transform、rank bitvectors 和稀疏 suffix-array samples 定位最长精确 seeds（Ferragina and Manzini, 2000）。当前实现使用 Ukkonen 在线 suffix-tree 算法生成 suffix order（Ukkonen, 1995）。方向一致性沿用 GeneMiner2 的 run-k 框架，其统计背景参见 Wald and Wolfowitz（1940）、Flajolet and Sedgewick（2009）及 Makri and Psillakis（2011）。可选的 alignment-shadow 模式仅对每个 locus 的有界 fragment 子集执行 seeded local alignment，并使用 affine gap costs 输出审计证据，不参与默认组装决策（Smith and Waterman, 1981；Gotoh, 1982）。

#### 2.3.2 可选的保守自动敏感招募

`--uce-recruit-mode auto` 是显式启用的两阶段策略；默认仍为 `fast`。第一阶段原样执行 `fast`，第二阶段只针对第一阶段没有选中 fragments 的 loci 再扫描一次 FASTQ，默认使用 recruitment *k*=21、step=1 和独立 verification *k*=19。未恢复 locus 子集只承担粗招募门控；对通过该门控的 fragments，程序再以完整 probe/reference 面板扩展候选并执行多 locus 歧义检查，避免因先缩小参考集合而把共享 reads 错判为 locus 特异证据。

敏感阶段要求 read pair 至少一端与目标 probe/reference 的局部比对覆盖不少于 45 bp 且 identity 不低于 80%，并只合并唯一支持一个未恢复 locus 的新增 FASTQ；快速阶段已经接受的 locus 和序列保持不变，合并后只执行一次组装。仅由敏感阶段恢复的 contig 还须满足长度至少 200 bp、目标 probe coverage 和 identity 均不低于 80%，且完整面板中不得存在满足相同 coverage/identity 条件且比对得分达到目标得分 95% 的其他 locus。该门控扩大的是候选 reads 的招募范围，不构成 orthology 或组装正确性的证明。

为保证可审计性，流程保留 `uce_filter_summary.fast.tsv`，在实际运行敏感阶段时另写 `uce_filter_summary.fallback.tsv`；`uce_recruit_passes.tsv` 记录每个 locus 的 fast、fallback 和最终来源，`uce_recruit_contig_probe_gate.tsv` 记录 fallback-only contig 的最终检查证据，被拒绝的结果归档至 `fallback_probe_rejected/` 而非静默删除。

#### 2.3.3 自适应逐 locus 选择

低证据 locus 默认不下采样。当合格候选少于 512 个、估计深度不高于 160×，或精确 seeds 覆盖少于 64 个参考区间中的 48 个时，reads 直接通过或回退至兼容选择器。对于饱和且广泛覆盖参考的 locus，UCEFilter 以约 80×为目标选择 core fragments，同时至少保留合格 fragments 的 60%。选择过程优先维持参考位置多样性；跨越参考末端的 fragments 按 overhang 长度分层，每个延伸长度保留少量最佳证据，从而保存由核心通向侧翼的 overlap ladder。

64-bin 位置覆盖、60% 最低保留比例、paired-fragment 原子选择和末端 overhang ladder 均为本文提出的 TStools 机制，不强行附会外部算法来源；其作用通过消融实验检验。

#### 2.3.4 参考锚定的加权图组装与质量控制

`uce-rust` 从每个 locus 的选择后 reads 构建参考锚定的加权 de Bruijn 图（Idury and Waterman, 1995；Xie et al., 2024；Yu et al., 2026）。节点权重综合 reads 深度与参考位置支持；paired fragments 为分支提供独立证据。默认 backbone 策略在每个 bubble 上进行一次有界前瞻并提交最佳分支，同时保留不依赖 paired-end 辅助的 core graph 作为安全候选。paired-fragment 分支支持、有界前瞻和 core-graph safeguard 为 TStools 专用机制。

候选 contig 依据 reads 支持而非单纯长度接受或拒绝。输出指标包括 unique-read density、slice-supported breadth、最大无支持 gap、k-mer depth coefficient of variation、maximum/median depth ratio、fragment support 和 paired-fragment support。【待补：正式分析采用的各 QC 阈值；若使用默认值，列入补充表 S1。】

#### 2.3.5 受证据约束的 rescue

可选 rescue 借鉴 MITObim 的 baiting-and-iterative-mapping 思路，也与 aTRAM 的 target-restricted assembly 范式相关（Hahn et al., 2013；Allen et al., 2015），但不复制参考序列以填补 gap。rescue 仅在 primary contig 已通过 QC 后运行，默认以固定 *k*=21 最多执行两轮：第一轮以原参考和已接受 contig 招募 reads；仍在增长的 locus 第二轮仅使用 contig 末端窗口。

每一侧新增序列必须至少为 30 bp，支持 breadth 不低于 85%，最大无支持 gap 不超过 30 bp，并具有至少两个独立 fragments，其中至少一个跨越冻结 core 与新增序列的边界。任一侧不满足条件即单侧回滚；若出现新的至少 150-bp 精确倒置重复，则整轮回滚。【待补：确认论文实验是否全部采用这些默认阈值。】rescue 的目标是增加有证据侧翼，而非最大化 contig 长度。

### 2.4 数据集与比较方法

#### 2.4.1 已知真值模拟

使用【待补：真实或模拟祖先参考、locus 数和长度分布】生成 UCE 核心及侧翼。序列沿【待补：树拓扑和枝长】用 Seq-Gen 模拟演化（Rambaut and Grass, 1997），再用 NGSNGS 生成 paired Illumina-like reads（Henriksen et al., 2023）。测试覆盖度设为 1×、3×、5×、10×、20×和 50×，参考—目标 divergence 设为 0.05、0.10、0.15和 0.20；每种组合重复【待补：重复数】次。另加入【待补：heterozygosity、indel、repeat、contamination 和 insert-size 设置】以评估困难场景。

比较 TStools、GeneMiner2、PHYLUCE + SPAdes、HybPiper 和 Captus + MEGAHIT。【待补：各软件版本与完整参数；保证线程数、输入和参考集合一致。】PHYLUCE、HybPiper 与 Captus 分别代表 assembly-first、target-bin 和 whole-sample assembly 路线（Johnson et al., 2016；Faircloth, 2016；Ortiz et al., 2026）。标准化组装统计可由 QUAST 提供（Gurevich et al., 2013），但核心准确性直接按 recovered locus 与已知目标序列的对应关系计算，包括 nucleotide identity、reference-covered fraction、完整正确 locus 比例、substitution/indel rate、chimera rate、错误延伸长度和未支持序列比例。

#### 2.4.2 v1.6.1 开发期回归

PR #12 对 `fast` 与 `auto` 进行了三组开发期回归。蜜蜂（*Apis mellifera*）和河豚（*Takifugu flavidus*）候选序列与同个体高质量基因组比对，用于判断新增结果是否由目标 locus 支持；这一设计提供强外部参照，但不是从模拟 reads 到目标序列的完整已知真值。另对两个珊瑚 WGS 样本 CRR2698937 和 CRR2698939 比较 fast/auto 的接受 locus 数、fast locus 丢失数和共有序列变化数。所有回归同时记录 FASTQ spill、wall time 和峰值 RSS。【待补：归档输入 accession、probe set、完整命令、同个体基因组版本、硬件和逐 locus 原始证据；正式投稿前复跑并固定校验和。】

#### 2.4.3 消融实验

为确定各机制的独立贡献，比较：（i）`fast` 与 `auto` 招募；（ii）legacy selection 与 adaptive selection；（iii）paired-fragment 原子保留开/关；（iv）位置分层选择开/关；（v）terminal overhang ladder 开/关；（vi）paired-end branch support 开/关；（vii）bounded look-ahead 开/关；（viii）rescue 开/关。主要响应变量为正确完整度、错误率、chimera rate、侧翼净增长、峰值 RSS 和 wall time。`auto` 另报告新增 locus 的 target-probe gate 通过率、panel near-tie 拒绝率以及相对 `fast` 的 locus 丢失与共有序列变化。【待补：统计模型、效应量和多重检验策略。】

#### 2.4.4 真实 UCE 数据与下游可用性

选用【待补：类群、样本数、BioProject/SRA accession、probe set 和测序类型】。为避免把真实数据中的树一致性误称为准确率，真实数据只报告 recovered loci、contig length、sample/locus occupancy、缺失率、variable sites、parsimony-informative sites 和 gene-tree/species-tree congruence。每个 locus 用 MAFFT 比对并以 trimAl 修剪（Katoh and Standley, 2013；Capella-Gutiérrez et al., 2009）；用 IQ-TREE 2 推断基因树，并在适用时用 ASTRAL-III 汇总物种树（Minh et al., 2020；Zhang et al., 2018）。【待补：实际命令、模型、bootstrap 和过滤阈值。】

#### 2.4.5 轻量化与跨平台基准

在三类环境运行同一完整 UCE 数据集：（i）Windows 11、【待补：CPU/内存】；（ii）macOS【待补：版本、Intel 或 Apple Silicon、CPU/内存】；（iii）Linux【待补：发行版、CPU/内存】，并对每个进程施加 16-GB 内存上限。记录完成率、wall time、峰值 RSS、CPU time、输入/输出量、最大临时磁盘和安装步骤数。与 PHYLUCE + SPAdes、HybPiper 和 Captus + MEGAHIT 的资源比较须使用同一 reads、参考、线程数及可比输出范围。

只有 Windows/macOS CI、发行二进制和实机基准全部通过后，正文才使用“支持普通 Windows 和 macOS 笔记本”的确定表述；此前保留为设计目标。

## 3 结果

### 3.1 核心 reads 招募保持兼容并降低资源消耗

在 1,000,000 对真实 UCE reads 和 3,579 个参考 loci 上，TStools Rust MainFilter 与历史 GeneMiner2 Haxe/C++ MainFilter 均产生 2,228 个非空 locus 文件，总输出均为 460,812,556 bytes；逐文件内容及逐 locus reads 计数一致。首次构建 cache 时，wall time 从 20.32 s 降至 3.52 s，峰值 RSS 从 381 MiB 降至 187 MiB；复用 cache 时，wall time 从 18.35 s 降至 3.72 s，峰值 RSS 从 375 MiB 降至 186 MiB。对应运行时间改善为 4.9–5.8 倍，内存降低约 50%。

【待补：硬件、操作系统、编译器、重复次数、均值或中位数及离散度。】  
【待补：完整 `refilter → original-rust` contig 与 GeneMiner2 基线的一致性及资源结果。】

| 情形 | 历史 GeneMiner2 MainFilter | TStools Rust MainFilter | 改善 |
| --- | ---: | ---: | ---: |
| 首次 cache + 过滤 | 20.32 s；381 MiB | 3.52 s；187 MiB | 5.8×；RSS −51.0% |
| 复用 cache + 过滤 | 18.35 s；375 MiB | 3.72 s；186 MiB | 4.9×；RSS −50.3% |

### 3.2 自动敏感招募的开发期回归

在同个体基因组外部参照检查中，`auto` 未丢失任何 `fast` 已接受 locus，也未改变 fast/auto 共有 locus 的序列。蜜蜂的接受 locus 数由 2,373 增至 2,419，新增 46 个；其中 45 个可由同个体基因组评分，44 个（97.78%）获得目标 locus 支持。河豚的接受 locus 数由 1,258 增至 1,281，新增 23 个；其中 21 个可评分，21 个（100%）获得目标 locus 支持。相对于参考 locus 总数的恢复比例分别由 91.77% 增至 93.50%，以及由 94.50% 增至 96.13%。两组均未观察到新增序列的倒置或 split alignment，但这些指标仍不是已知真值下的 substitution、indel 或 chimera rate。

| 数据集 | Fast accepted | Auto accepted | 新增 | 新增可评分 | 新增目标 locus 支持 | Fast 丢失 / 共有序列改变 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| *Apis mellifera* | 2,373 | 2,419 | 46 | 45 | 44（97.78%） | 0 / 0 |
| *Takifugu flavidus* | 1,258 | 1,281 | 23 | 21 | 21（100%） | 0 / 0 |

在珊瑚 WGS 回归中，CRR2698937 的接受 locus 数由 719 增至 892，CRR2698939 由 1,242 增至 1,855；分别新增 173 和 613 个 locus，且均无 fast locus 丢失或共有序列改变。两样本 `auto` 运行的 wall time 为 10 min 56 s，峰值 RSS 为 261 MiB，四次 FASTQ 扫描均报告 0.0 MiB spill。【待补：运行硬件、输入规模、版本锁定和逐 locus 外部验证；因此该资源结果只作为开发期工程回归，不与其他流程作正式性能比较。】

| 珊瑚样本 | Fast accepted | Auto accepted | 新增 | Fast 丢失 | 共有序列改变 |
| --- | ---: | ---: | ---: | ---: | ---: |
| CRR2698937 | 719 | 892 | 173 | 0 | 0 |
| CRR2698939 | 1,242 | 1,855 | 613 | 0 | 0 |

### 3.3 UCE 恢复的准确性与完整度

在已知真值模拟中，TStools 在【待补：覆盖度与 divergence 范围】获得【待补：median identity】的 nucleotide identity、【待补】的正确完整 locus 比例和【待补】的 chimera rate。相较 PHYLUCE + SPAdes、HybPiper、Captus + MEGAHIT 和 GeneMiner2，TStools 的主要差异为【待补：效应量与置信区间】，而不是单纯恢复更多 loci。

【待补：低深度结果。】  
【待补：高 divergence 结果。】  
【待补：重复、杂合和污染场景结果。】

### 3.4 自适应选择与 rescue 的独立贡献

开启 adaptive selection 后，饱和 locus 的输入 fragments 从【待补】降至【待补】，峰值 RSS 和 wall time分别降低【待补】与【待补】，而正确完整度变化为【待补】。关闭位置分层或 terminal ladder 后，【待补：侧翼或错误率变化】，说明【待补：机制解释】。

rescue 使通过 QC 的 contig 中位长度由【待补】增加至【待补】，新增侧翼的【待补】%具有完整 reads 支持；substitution/indel 和 chimera rate 分别变化【待补】。被 side-specific rollback 或 inverted-repeat guard 拒绝的 loci 为【待补】个。

### 3.5 普通笔记本可完成 UCE 分析

在 16-GB Windows 笔记本上，TStools 完成【待补：数据规模】需要【待补：时间】，峰值 RSS 为【待补】，临时磁盘为【待补】；在 macOS 笔记本上对应结果为【待补】。两种平台均【待补：是否无需 Python、SPAdes 和服务器级硬件完成全流程】。与其他流程相比，TStools 的资源优势为【待补：定量结果】，其序列准确性和完整度差异见第 3.3 节。

### 3.6 真实数据的下游可用性

TStools 在【待补：样本数】中恢复中位【待补】个 loci，矩阵 occupancy 为【待补】，每 locus 的 parsimony-informative sites 为【待补】。与比较流程相比，其 gene-tree/species-tree congruence 为【待补】。这些指标只反映真实数据中的恢复和下游一致性，不作为已知真值准确率。

## 4 讨论

### 4.1 从“更多 contigs”转向“单位计算成本下的可靠序列”

TStools 的核心设计不是增加一组彼此独立的功能，而是缩小问题：先用目标参考招募可能相关的 paired fragments，再以位置、方向和精确匹配证据控制每个 locus 的计算规模，最后只接受 reads 支持的图路径和侧翼延伸。默认 fast 路径维持一次扫描；当完整跨物种参考的分化侧翼使 *k*=31 招募不足时，`auto` 只扩大未恢复 loci 的候选范围，并在完整面板、read 局部比对和 fallback-only contig 三层门控后接纳新增结果。这一设计把 UCE 恢复从通用组装问题转化为 target-specific、evidence-bounded 的恢复问题。

Bossert et al.（2024）说明，恢复 locus 数量最多的流程未必同时得到缺失最少、信息位点最多或 gene-tree congruence 最高的矩阵。因此，TStools 的主要评价标准应是 known-truth accuracy、usable sequence 和 resource cost 的联合表现，而不是以 contig 数量或长度作为唯一终点。Captus 的全样本 MEGAHIT 路线可充分利用 off-target reads，且其组装成本不随 target 数量线性增加；HybPiper 的逐 target 路线则便于基因级提取和旁系同源检查（Johnson et al., 2016；Ortiz et al., 2026）。TStools 与两者的差异在于更早限制计算范围并将证据边界直接纳入 UCE 组装，因此应被定位为不同数据规模和硬件条件下的补充选择，而非宣称普遍替代所有 assembly-first 流程。

### 4.2 与 GeneMiner2 的继承关系和新增贡献

TStools 继承 GeneMiner/GeneMiner2 的参考引导过滤、run-based 精细选择和加权 de Bruijn 图思想（Xie et al., 2024；Yu et al., 2026）。这些来源应在摘要、方法和讨论中明确引用。本文新增贡献应限定为 UCEFilter 默认单次扫描架构、paired-fragment 原子存储与选择、64-bin 位置覆盖控制、terminal overhang ladder、paired-end branch support、bounded look-ahead、core-graph safeguard、逐侧 rollback 的 evidence-only rescue，以及 v1.6.1 中只针对未恢复 loci 的保守敏感招募、完整面板歧义检查和 fallback-only contig gate。每项新增机制均需由消融实验而非描述性语言证明价值。

### 4.3 轻量化是可复现实验结论，而不是对 SPAdes 的泛化判断

TStools 不为每个样本或 locus 调用通用组装器，而只对自适应选择后的 locus-specific reads 建图，因此有望降低峰值内存、临时存储和安装依赖。论文不应写成“SPAdes 必须使用超大内存服务器”；SPAdes 的实际资源需求取决于数据规模和参数。更准确的表述是：**TStools 的核心 UCE 路径不依赖 Python 或 SPAdes，并被设计为在普通个人电脑上完成常见规模的数据分析。**完成跨平台实测后，可进一步写为：**完整 UCE 分析在普通 16-GB Windows 和 macOS 笔记本上完成，无需服务器级硬件。**

### 4.4 局限性

首先，参考引导方法可能漏掉与 bait 高度分化的 reads，并可能在重复、旁系同源或高杂合区域产生歧义。`auto` 可恢复部分低敏感度 fast pass 漏失的 loci，但其完整面板歧义检查只覆盖用户提供的 probe/reference 集合，不能排除面板外旁系同源、重复或污染；read/contig gate 也不能证明所有输出均为正确 ortholog。完整面板 near-tie 检查还可能在 fallback contig 较多时主导运行时间。其次，短 reads 无法可靠解析超过 insert size 的完全重复或复杂结构；受证据约束的选择、rescue 和保守回滚会以牺牲部分敏感性换取较低的错误风险。第三，同个体基因组比对和真实数据的 gene-tree congruence 均不是已知真值准确率，最终结论必须以已知真值模拟和充分的类群采样支持。

最后，改进计算恢复并不能消除 UCE marker 本身的 ascertainment bias。UCE loci 在基因组中并非随机分布，在存在广泛 introgression 或 post-speciation gene flow 时，UCE 数据可能偏离物种树并压缩枝长（Foley and Murphy, 2026）。因此，TStools 解决的是 reads 到序列的计算恢复问题，而不是 UCE 位点选择、重组、基因流或物种树模型的生物学局限。

## 5 结论

TStools 以默认单次扫描的 UCE 证据选择、仅面向未恢复 loci 的可选保守敏感扫描、参考锚定的加权图组装和只接受 reads 支持的 rescue，减少了 UCE 恢复对通用组装器和高内存计算环境的依赖。v1.6.1 的开发期回归表明，`auto` 可在不丢失 fast 已接受 loci、也不改变共有序列的情况下增加候选恢复，但该结果仍需已知真值和正式跨流程实验确认。【待补：在正式结果支持后加入一句定量结论。】其价值不在于无条件恢复最长或最多的 contigs，而在于以较少计算资源获得正确、可审计且适合下游分析的 UCE 序列。

## 数据和代码可用性

- TStools 源代码：https://github.com/GUIBA-EX/TStools。
- 本文软件快照：v1.6.1，Git commit `3d198e3`，https://github.com/GUIBA-EX/TStools/releases/tag/v1.6.1；【待补：Zenodo/其他永久归档 DOI】。
- 自动敏感招募的实现与开发期验证记录：https://github.com/GUIBA-EX/TStools/pull/12。
- 模拟脚本、参数及 truth sequences：【待补：地址】。
- 真实数据 accession：【待补：BioProject/SRA】。
- 完整 benchmark 输出与环境信息：【待补：地址】。

## 致谢、经费和利益冲突

**致谢：**【待补】  
**经费：**【待补】  
**作者贡献：**【待补：CRediT】  
**利益冲突：**【待补】

## 图表建议

- **图 1：** 仅画两条主流程：核心 `MainFilter → refilter → original-rust`；UCE `UCEFilter fast [→ unresolved-only auto fallback] → uce-rust → optional rescue`。突出 fast 默认一次扫描、auto 仅对未恢复 loci 再扫描、paired-fragment 原子保留和 evidence-only rescue。
- **图 2：** 覆盖度 × 参考 divergence 下，各工具的正确完整 locus 比例、identity 和 chimera rate。
- **图 3：** fast/auto、adaptive selection 与 rescue 消融：准确性、侧翼增益、峰值 RSS 和时间。
- **图 4：** Windows、macOS 和 Linux 的峰值 RSS、时间及临时磁盘。
- **表 1：** 数据集、软件版本、参数与硬件。
- **表 2：** known-truth 准确性、fast/auto 新增与拒绝证据、真实数据可用性和资源指标；不要只列 recovered loci。

## 参考文献

Allen JM, Huang DI, Cronk QC, Johnson KP. 2015. aTRAM—automated target restricted assembly method: a fast method for assembling loci across divergent taxa from next-generation sequencing data. *BMC Bioinformatics* 16:98. https://doi.org/10.1186/s12859-015-0515-2

Andermann T, Cano Á, Zizka A, Bacon C, Antonelli A. 2018. SECAPR—a bioinformatics pipeline for the rapid and user-friendly processing of targeted enriched Illumina sequences, from raw reads to alignments. *PeerJ* 6:e5175. https://doi.org/10.7717/peerj.5175

Bankevich A, Nurk S, Antipov D, et al. 2012. SPAdes: A new genome assembly algorithm and its applications to single-cell sequencing. *Journal of Computational Biology* 19:455–477. https://doi.org/10.1089/cmb.2012.0021

Bloom BH. 1970. Space/time trade-offs in hash coding with allowable errors. *Communications of the ACM* 13:422–426. https://doi.org/10.1145/362686.362692

Bossert S, Pauly A, Danforth BN, Orr MC, Murray EA. 2024. Lessons from assembling UCEs: A comparison of common methods and the case of *Clavinomia* (Halictidae). *Molecular Ecology Resources* 24:e13925. https://doi.org/10.1111/1755-0998.13925

Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics* 25:1972–1973. https://doi.org/10.1093/bioinformatics/btp348

Faircloth BC. 2016. PHYLUCE is a software package for the analysis of conserved genomic loci. *Bioinformatics* 32:786–788. https://doi.org/10.1093/bioinformatics/btv646

Faircloth BC, McCormack JE, Crawford NG, Harvey MG, Brumfield RT, Glenn TC. 2012. Ultraconserved elements anchor thousands of genetic markers spanning multiple evolutionary timescales. *Systematic Biology* 61:717–726. https://doi.org/10.1093/sysbio/sys004

Ferragina P, Manzini G. 2000. Opportunistic data structures with applications. In: *Proceedings of the 41st Annual Symposium on Foundations of Computer Science*. pp. 390–398. https://doi.org/10.1109/SFCS.2000.892127

Flajolet P, Sedgewick R. 2009. *Analytic Combinatorics*. Cambridge University Press. https://doi.org/10.1017/CBO9780511801655

Foley NM, Murphy WJ. 2026. Phylogenomic blind spots: the limits of UCE and BUSCO loci in the presence of gene flow. *Molecular Biology and Evolution* 43:msag155. https://doi.org/10.1093/molbev/msag155

Gotoh O. 1982. An improved algorithm for matching biological sequences. *Journal of Molecular Biology* 162:705–708. https://doi.org/10.1016/0022-2836(82)90398-9

Gurevich A, Saveliev V, Vyahhi N, Tesler G. 2013. QUAST: quality assessment tool for genome assemblies. *Bioinformatics* 29:1072–1075. https://doi.org/10.1093/bioinformatics/btt086

Hahn C, Bachmann L, Chevreux B. 2013. Reconstructing mitochondrial genomes directly from genomic next-generation sequencing reads—a baiting and iterative mapping approach. *Nucleic Acids Research* 41:e129. https://doi.org/10.1093/nar/gkt371

Henriksen RA, Zhao L, Korneliussen TS. 2023. NGSNGS: next-generation simulator for next-generation sequencing data. *Bioinformatics* 39:btad041. https://doi.org/10.1093/bioinformatics/btad041

Idury RM, Waterman MS. 1995. A new algorithm for DNA sequence assembly. *Journal of Computational Biology* 2:291–306. https://doi.org/10.1089/cmb.1995.2.291

Johnson MG, Gardner EM, Liu Y, et al. 2016. HybPiper: extracting coding sequence and introns for phylogenetics from high-throughput sequencing reads using target enrichment. *Applications in Plant Sciences* 4:1600016. https://doi.org/10.3732/apps.1600016

Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30:772–780. https://doi.org/10.1093/molbev/mst010

Knyshov A, Gordon ERL, Weirauch C. 2021. New alignment-based sequence extraction software (ALiBaSeq) and its utility for deep level phylogenetics. *PeerJ* 9:e11019. https://doi.org/10.7717/peerj.11019

Li D, Liu CM, Luo R, Sadakane K, Lam TW. 2015. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics* 31:1674–1676. https://doi.org/10.1093/bioinformatics/btv033

Makri FS, Psillakis ZM. 2011. On success runs of a fixed length in Bernoulli sequences: exact and asymptotic results. *Computers & Mathematics with Applications* 61:761–772. https://doi.org/10.1016/j.camwa.2010.12.023

Minh BQ, Schmidt HA, Chernomor O, et al. 2020. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution* 37:1530–1534. https://doi.org/10.1093/molbev/msaa015

Ortiz EM, Höwener A, Shigita G, et al. 2026. A novel phylogenomics pipeline reveals extensive topological conflict in the evolution of the angiosperm order Cucurbitales. *Systematic Biology*:syag046. https://doi.org/10.1093/sysbio/syag046

Putze F, Sanders P, Singler J. 2009. Cache-, hash-, and space-efficient Bloom filters. *ACM Journal of Experimental Algorithmics* 14. https://doi.org/10.1145/1498698.1594230

Rambaut A, Grass NC. 1997. Seq-Gen: an application for the Monte Carlo simulation of DNA sequence evolution along phylogenetic trees. *Bioinformatics* 13:235–238. https://doi.org/10.1093/bioinformatics/13.3.235

Raza M, Ortiz EM, Schwung L, Shigita G, Schaefer H. 2023. Resolving the phylogeny of *Thladiantha* (Cucurbitaceae) with three different target capture pipelines. *BMC Ecology and Evolution* 23:75. https://doi.org/10.1186/s12862-023-02185-z

Ribeiro C, Oliveira L, Batista R, De Sousa M. 2021. UCEasy: a software package for automating and simplifying the analysis of ultraconserved elements (UCEs). *Biodiversity Data Journal* 9:e78132. https://doi.org/10.3897/BDJ.9.e78132

Smith BT, Harvey MG, Faircloth BC, Glenn TC, Brumfield RT. 2014. Target capture and massively parallel sequencing of ultraconserved elements for comparative studies at shallow evolutionary time scales. *Systematic Biology* 63:83–95. https://doi.org/10.1093/sysbio/syt061

Smith TF, Waterman MS. 1981. Identification of common molecular subsequences. *Journal of Molecular Biology* 147:195–197. https://doi.org/10.1016/0022-2836(81)90087-5

Ukkonen E. 1995. On-line construction of suffix trees. *Algorithmica* 14:249–260. https://doi.org/10.1007/BF01206331

Wald A, Wolfowitz J. 1940. On a test whether two samples are from the same population. *The Annals of Mathematical Statistics* 11:147–162. https://doi.org/10.1214/aoms/1177731909

Xie P, Guo Y, Teng Y, Zhou W, Yu Y. 2024. GeneMiner: a tool for extracting phylogenetic markers from next-generation sequencing data. *Molecular Ecology Resources* 24:e13924. https://doi.org/10.1111/1755-0998.13924

Yu X, Tang Z, Zhang Z, Song Y, He H, Shi Y, Hou J, Yu Y. 2026. GeneMiner2: accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26:e70111. https://doi.org/10.1111/1755-0998.70111

Zhang C, Rabiee M, Sayyari E, Mirarab S. 2018. ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. *BMC Bioinformatics* 19(Suppl 6):153. https://doi.org/10.1186/s12859-018-2129-y

Zhang Z, Xie P, Guo Y, Zhou W, Liu E, Yu Y. 2022. Easy353: a tool to get Angiosperms353 genes for phylogenomic research. *Molecular Biology and Evolution* 39:msac261. https://doi.org/10.1093/molbev/msac261
