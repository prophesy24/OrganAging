# R包配置 ----
library(ggVennDiagram)
library(ggplot2)
library(readxl)
library(VennDiagram)
library(dplyr)
library(tidyr)
library(ggalluvial)  # 桑基图/冲击图绘制
library(networkD3)
library(readr)
library(tibble)  # 提供column_to_rownames()
library(pheatmap)
library(ggraph)
library(tidygraph)
library(STRINGdb)
library(igraph)
library(tidyverse)
library(shinybody)
library(scales)
library(htmlwidgets)
library(circlize)
library(RColorBrewer)
library(ggalluvial)
library(cowplot)

# 将SXF血浆蛋白与marker做交集
BIOMARKER <- read_excel("./01_data/BIOMARKER.xlsx")
View(BIOMARKER)
SXF_protein <- read_excel("./01_data/SXF protein.xlsx")
View(SXF_protein)
BIOMARKER <- as.data.frame(BIOMARKER)
SXF_protein <- as.data.frame(SXF_protein)

# 提取基因列表并去重
organ_aging <- unique(BIOMARKER$Gene)
SXF_protein <- unique(SXF_protein$Gene)

# 获取共有基因
common_genes <- intersect(organ_aging, SXF_protein)
print(common_genes)

common_genes <- as.data.frame(common_genes)

# Venn图 

# 获取两个集合
set1 <- organ_aging
set2 <- SXF_protein

# 计算交集
intersect_len <- length(intersect(set1, set2))

# 创建 Venn 图
venn.plot <- draw.pairwise.venn(
  area1 = length(set1),
  area2 = length(set2),
  cross.area = intersect_len,
  category = c("Organ aging biomarker", "Plasma protein"),
  fill = c("#657DAF", "#A7687D"),
  alpha = 0.6,
  cat.pos = c(0, 0), # 类别标签的位置
  cat.dist = c(0.03, 0.03),
  cat.cex = 1.2,
  cex = 1.5,
  fontface = "bold",
  lty = "blank"
)

# 保存图像到文件
pdf("Venn_organ_vs_plasma.pdf", width = 6, height = 6)
grid.draw(venn.plot)
dev.off()

# 输出共有基因
write.csv(data.frame(Common_Genes = common_genes), 
          "common_genes.csv", 
          row.names = FALSE)

# 与HPA的分泌蛋白overlap ----

library(hpar)
# 获取HPA所有蛋白质分类数据（含分泌蛋白标记）
secretome <- read.table("./01_data/sa_location_Secreted.tsv", header=TRUE, sep="\t")


# 取交集：筛选你的标志物中属于分泌蛋白的基因
overlap_genes <- inner_join(BIOMARKER, secretome, by = "Gene")

# 保存结果
write.csv(overlap_genes, "secretome_biomarkers_overlap.csv", row.names = FALSE)

# 统计结果
cat(sprintf("你的标志物中有%d个分泌蛋白（占比%.1f%%）",
            nrow(overlap_genes),
            nrow(overlap_genes)/nrow(BIOMARKER)*100))


common_genes$group <- "plasma aging marker"
overlap_genes$group <- "secreted aging marker"


secreted <- secretome$Gene


# combined
# 将全部overlap
# combined
venn_genes <- unique(c(organ_aging, SXF_protein, secreted))
venn.plot <- draw.triple.venn(
  area1 = length(organ_aging),
  area2 = length(SXF_protein),
  area3 = length(secreted),
  n12 = length(intersect(organ_aging, SXF_protein)),
  n23 = length(intersect(SXF_protein, secreted)),
  n13 = length(intersect(organ_aging, secreted)),
  n123 = length(Reduce(intersect, list(organ_aging, SXF_protein, secreted))),
  category = c("Organ Aging", "Plasma Protein", "Secreted"),
  fill = c("#FC8D62", "#8DA0CB", "#66C2A5"),
  alpha = 0.6,
  lty = "blank",
  cat.cex = 1.2,
  cex = 1.5
)
# 保存为 PDF
pdf("Triple_Venn.pdf", width = 7, height = 7)
grid.draw(venn.plot)
dev.off()

All_overlap_genes <- Reduce(intersect, list(organ_aging, SXF_protein, secreted))
overlap_aging_Plasma_secreted <- All_overlap_genes
# res output
overlap_aging_Plasma_secreted <- data.frame(overlap_aging_Plasma_secreted)
write.xlsx(overlap_aging_Plasma_secreted, file = "./result/All overlap.xlsx")

#将基因映射到器官 ----
# 安装并加载hpar包

# 获取HPA正常组织表达数据（含分泌蛋白信息）
# 使用getHpa()函数获取分泌蛋白数据
#组织数据
hpa_tissue <- read.delim("./01_data/normal_ihc_data.tsv")

# 分泌蛋白数据
hpa_secretome <- read.delim("./01_data/proteinatlas.tsv") 
# 步骤2：读取数据
hpa_secretome <- secretome
# 验证列名
colnames(hpa_tissue)      # 应包含 Gene, Tissue, Cell.type, Level
colnames(hpa_secretome)   # 应包含 Gene, Secretory.score

# 统计不同分泌定位的蛋白数量
location_count <- hpa_secretome %>%
  filter(!is.na(Secretome.location)) %>%
  count(Secretome.location) %>%
  arrange(desc(n))
location_count$percent <- round(location_count$n / sum(location_count$n) * 100, 1)
location_count$label <- paste0(location_count$Secretome.location, "\n", location_count$percent, "%")
# HPA数据中分泌蛋白的定位
ggplot(location_count, aes(x = "", y = n, fill = Secretome.location)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  labs(title = "分泌蛋白的亚定位分布（饼图）", fill = "分泌定位") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "Set3")

ggsave("HPA_Secretome_location_pie.pdf", width = 6, height = 6, dpi = 300)



# 将我的数据映射到器官 ----

# 从HPA加载组织表达数据（若已有hpa_tissue数据框则跳过）

# 跨器官衰老网络分析 ----
# 安装和加载必要包



# Step 1: 读取你的基因列表（Gene_symbol + organ）
gene_list <- read_excel("./01_data/The cross-organ aging network.xlsx")
View(gene_list)

# Step 2: 加载 protein.aliases 文件，手动映射 Gene_symbol -> STRING ID
aliases <- read.delim("./01_data/9606.protein.aliases.v11.5.txt.gz", header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

# 过滤常用命名类型，选取最有代表性的映射
aliases_filtered <- aliases %>% filter(source %in% c("Ensembl_HGNC", "Ensembl", "Gene_Name"))
gene_map <- merge(gene_list, aliases_filtered, by.x="Gene_symbol", by.y="alias")
colnames(gene_map)[colnames(gene_map) == "stringId"] <- "STRING_id"

# Step 3: 加载 protein.links 文件，提取我们需要的 PPI 关系
links <- read.delim("./01_data/9606.protein.links.v11.5.txt.gz", header=TRUE, sep=" ")

# 只保留我们的输入基因之间的交互
input_ids <- gene_map$X.string_protein_id
ppi_filtered <- links %>%
  filter(protein1 %in% input_ids & protein2 %in% input_ids)

# Step 4: 给交互表加上器官信息
ppi_annotated <- ppi_filtered %>%
  left_join(gene_map[, c("X.string_protein_id", "Organ")], by = c("protein1" = "X.string_protein_id")) %>%
  rename(organ_from = Organ) %>%
  left_join(gene_map[, c("X.string_protein_id", "Organ")], by = c("protein2" = "X.string_protein_id")) %>%
  rename(organ_to = Organ)

# Step 5: 提取跨器官的交互
ppi_cross_organ <- ppi_annotated %>% filter(organ_from != organ_to)

# Step 6: 构建器官-器官连接表
organ_edges <- ppi_cross_organ %>%
  group_by(organ_from, organ_to) %>%
  summarise(interaction_count = n()) %>%
  ungroup()

# Step 7: 构建器官网络图
organ_graph <- graph_from_data_frame(organ_edges, directed = FALSE)

# Step 8: 可视化
plot(
  organ_graph,
  vertex.size = 30,
  vertex.label.cex = 1.2,
  edge.width = E(organ_graph)$interaction_count / max(E(organ_graph)$interaction_count) * 10,
  main = "Organ-to-Organ Interaction Network via PPI"
)

# 计算每个节点的 strength（即加权度）
strengths <- strength(organ_graph, weights = E(organ_graph)$interaction_count)

# 对 strength 进行缩放，映射为合适的节点大小（比如 10~40）
scaled_strength <- scales::rescale(strengths, to = c(10, 40))

set.seed(123)
plot(
  organ_graph,
  layout = layout_with_fr(organ_graph),  # Fruchterman-Reingold 布局
  vertex.size = scaled_strength,
  vertex.label.cex = 1.2,
  vertex.label.color = "black",
  vertex.color = "skyblue",
  edge.width = scales::rescale(E(organ_graph)$interaction_count, to = c(1, 10)),
  edge.color = "gray70",
  main = "Organ PPI Network (Node size = Strength)"
)

# 中心性分析 ----
# Degree (连接了多少其他器官)
deg <- igraph::degree(organ_graph, mode = "all")


# Strength (总边权和 = 所有交互蛋白数)
strength <- strength(organ_graph, mode = "all", weights = E(organ_graph)$interaction_count)

# Betweenness
btw <- betweenness(organ_graph, weights = 1/E(organ_graph)$interaction_count, normalized = TRUE)

# Closeness
cls <- closeness(organ_graph, weights = 1/E(organ_graph)$interaction_count, normalized = TRUE)

# Eigenvector
eig <- eigen_centrality(organ_graph, weights = E(organ_graph)$interaction_count)$vector

# 组合成表格
organ_centrality <- data.frame(
  Organ = names(deg),
  Degree = deg,
  Strength = strength,
  Betweenness = btw,
  Closeness = cls,
  Eigenvector = eig
)

# 查看 Top hub organs
organ_centrality %>%
  arrange(desc(Strength)) %>%
  head(10)

write.csv(organ_centrality, file = "./Hub organ rank.csv")



# 再次绘图

# 添加中心性信息
organ_graph <- set_vertex_attr(organ_graph, "Strength", value = strength)

# 转换为 tidygraph 对象
tg <- as_tbl_graph(organ_graph)
# 创建颜色映射
strength_range <- range(tg %>% activate(nodes) %>% pull(Strength))
tg <- tg %>% 
  mutate(
    strength_scaled = (Strength - strength_range[1]) / diff(strength_range),
    label = name
  )

# 开始绘图
ggraph(tg, layout = "fr") + 
  geom_edge_link(aes(width = interaction_count), color = "gray80", alpha = 0.7) +
  geom_node_point(aes(size = Strength, color = strength_scaled)) +
  geom_node_text(aes(label = label), repel = TRUE, size = 4, color = "black") +
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  scale_size(range = c(5, 15)) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  ggtitle("Organ-level PPI Network of Aging Markers")

# 只显示重要节点
# 只对 Strength 排名前10 的节点加 label
tg <- tg %>%
  mutate(
    show_label = ifelse(rank(-Strength) <= 10, name, "")
  )

ggraph(tg, layout = "fr") + 
  geom_edge_link(aes(width = interaction_count), color = "gray80", alpha = 0.7) +
  geom_node_point(aes(size = Strength, color = strength_scaled)) +
  geom_node_text(aes(label = show_label), repel = TRUE, size = 4, color = "black") +
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  scale_size(range = c(5, 15)) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  ggtitle("Organ-level PPI Network of Aging Markers")


# 器官映射图片交互网络 ----
# 加载内置的男性人体器官坐标（human male）
# 加载内置解剖图 key 数据（人体器官坐标）
# 加载数据
library(gganatogram)
data("hgMale_key", package="gganatogram") # 加载解剖模板&#8203;:contentReference[oaicite:1]{index=1}

organ_centrality <- organ_centrality %>%
  mutate(strength_scaled = rescale(Strength))

# 器官映射表（之后要修改）
organ_name_map <- tibble::tribble(
  ~Organ,           ~organ,
  "Adipose tissue", "adipose_tissue",
  "Bladder",        "urinary_bladder",
  "Bone",           "bone",
  "Bone marrow",    "bone_marrow",
  "Brain",          "brain",
  "Breast",         "breast",
  "Gingiva",        "tongue",  # 近似代表
  "Heart",          "heart",
  "Intestinal",     "small_intestine",  # 可选 "colon" / "duodenum"
  "Kidney",         "kidney",
  "Liver",          "liver",
  "Lung",           "lung",
  "Muscle",         "skeletal_muscle",  # or smooth_muscle
  "Ovary",          NA,  # hgMale_key 是 male，没有 ovary
  "PBMC",           "leukocyte",
  "Pancreatic",     "pancreas",
  "Plasma",         NA,  # 无法映射
  "Skin",           "skin",
  "Spleen",         "spleen",
  "Testis",         "testis",
  "Vascular",       "aorta"  # 代表大血管
)


library(gganatogram)
library(dplyr)
library(ggplot2)

# 添加 organ 名称列用于匹配
organ_centrality_mapped <- organ_centrality %>%
  left_join(organ_name_map, by = "Organ") %>%
  filter(!is.na(organ)) %>%
  mutate(strength_scaled = scales::rescale(Strength)) %>%
  select(organ, strength_scaled)

# 合并到 hgMale_key
anat_data <- hgMale_key %>%
  select(-value) %>%  # 如果存在 value 列，先移除
  left_join(organ_centrality_mapped, by = "organ") %>%
  rename(value = strength_scaled)


# 绘图
gganatogram(data = anat_data,
            fill = "value",
            organism = "human",
            sex = "male",
            outline = TRUE,
            fillOutline = "gray90") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "YlOrRd"),
                       na.value = "white") +
  theme_void() +
  ggtitle("Organ Strength (Hubness) in Cross-Organ Aging Network")


# shinybody绘图 ----
# 安装并加载shinybody包（假设已发布在CRAN或GitHub）
# devtools::install_github("shinybody/shinybody")  # 如果从GitHub安装
# 使用示例数据进绘图
organ_data_male <- data.frame(
  organ = c("brain", "testis", "liver", "lung", "kidney", "pancreas", "adipose_tissue",
            "skeletal_muscle", "bone", "bone_marrow", "heart", "colon",
            "tongue", "bladder", "skin", "aorta", "breast", "spleen", "circulatory_system"),
  value = c(0.58, 0.44, 0.40, 0.07, 0.63, 0.34, 0.31, 0.31, 0.24, 0.19, 0.17, 0.14, 0.13, 0.13, 0.13, 0.10, 0.10, 0.01, 1),
  hovertext = c("Brain", "Testis", "Liver", "Lung", "Kidney", "Pancreatic",
                "Adipose_tissue","Muscle", "Bone", "Bone_marrow", "Heart", "Intestine",
                "Gingiva", "Bladder", "Skin", "Artery", "Breast", "Spleen", "Plasma")
)



# 为每个器官分配颜色
organ_data_male$color <- scales::col_numeric(
  palette = "YlOrRd",
  domain = range(organ_data_male$value)
)(organ_data_male$value)
# 使用scales::col_numeric来根据value字段进行颜色分配

# 创建交互式人体图
human(
  gender = "male",  # 或 "female"
  organ_df = organ_data_male,
  select_color = "yellow",
  width = 800,
  height = 1000
)

# female organ PPI
organ_data_female <- data.frame(
  organ = c("brain", "ovary", "liver", "lung", "kidney", "pancreas", "adipose_tissue",
            "skeletal_muscle", "bone", "bone_marrow", "heart", "colon",
            "tongue", "bladder", "skin", "aorta", "breast", "spleen", "circulatory_system"),
  value = c(0.58, 0.29, 0.40, 0.07, 0.63, 0.34, 0.31, 0.31, 0.24, 0.19, 0.17, 0.14, 0.13, 0.13, 0.13, 0.10, 0.10, 0.01, 1),
  hovertext = c("Brain", "Ovary", "Liver", "Lung", "Kidney", "Pancreatic",
                "Adipose_tissue","Muscle", "Bone", "Bone_marrow", "Heart", "Intestine",
                "Gingiva", "Bladder", "Skin", "Artery", "Breast", "Spleen", "Plasma")
)


# 为每个器官分配颜色
organ_data_female$color <- scales::col_numeric(
  palette = "YlOrRd",
  domain = range(organ_data_female$value)
)(organ_data_female$value)
# 使用scales::col_numeric来根据value字段进行颜色分配

# 创建交互式人体图
human(
  gender = "female",  # 或 "female"
  organ_df = organ_data_female,
  select_color = "yellow",
  width = 800,
  height = 1000
)

# 将图像输出

# 保存男性图
p_male <- human(
  gender = "male",
  organ_df = organ_data_male,
  select_color = "yellow",
  width = 800,
  height = 1000
)
saveWidget(p_male, "male_ppi.html", selfcontained = TRUE)

# 保存女性图
p_female <- human(
  gender = "female",
  organ_df = organ_data_female,
  select_color = "yellow",
  width = 800,
  height = 1000
)
saveWidget(p_female, "female_ppi.html", selfcontained = TRUE)


--------------

# 自定义器官颜色映射函数

map_strength_to_color <- function(strength) {
  colors <- colorRampPalette(c("#FFF7EC", "#FDBB84", "#EF6548"))(100)
  colors[cut(strength, breaks = seq(0, 1, length.out = 101))]
}

# 器官映射表（根据shinybody的器官ID调整）
organ_name_map <- tibble::tribble(
  ~Organ,           ~organ_id,
  "Adipose tissue", "adipose",
  "Bladder",        "bladder",
  "Bone",           "bone",
  "Bone marrow",    "bone_marrow",
  "Brain",          "brain",
  "Breast",         "breast",
  "Gingiva",        "tongue",
  "Heart",          "heart",
  "Intestinal",     "intestine",
  "Kidney",         "kidney",
  "Liver",          "liver",
  "Lung",           "lung",
  "Muscle",         "muscle",
  "Ovary",          "ovary",
  "PBMC",           "blood_cells",
  "Pancreatic",     "pancreas",
  "Plasma",         "blood",
  "Skin",           "skin",
  "Spleen",         "spleen",
  "Testis",         "testis",
  "Vascular",       "artery"
)

# 数据预处理
organ_centrality_mapped <- organ_centrality %>%
  left_join(organ_name_map, by = "Organ") %>%
  filter(!is.na(organ_id)) %>%
  mutate(
    strength_scaled = scales::rescale(Strength),
    color = map_strength_to_color(strength_scaled)
  )

# 绘制男性人体图
male_plot <- human(
  gender = "male",
  organ_df = organ_centrality_mapped,
  select_color = "color",  # 使用颜色填充
  width = 600,
  height = 800   # 显示颜色图例
)

# 绘制女性人体图（注意卵巢等性别特异性器官）
female_plot <- human(
  gender = "female",
  organ_df = organ_centrality_mapped %>%
    mutate(organ_id = ifelse(organ_id == "testis", "ovary", organ_id)),
  select_color = "color",  # 使用颜色填充
  width = 600,
  height = 800# 显示颜色图例
)

# 双图并列输出
library(patchwork)
combined_plot <- male_plot + female_plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

print(combined_plot)

#看器官对之间的所有交互蛋白 ----
# 设置器官对
organA <- "Heart"
organB <- "Liver"

# 提取两者之间的交互（注意双向）
organ_pair_ppi <- ppi_cross_organ %>%
  filter((organ_from == organA & organ_to == organB) |
           (organ_from == organB & organ_to == organA))

# 加回蛋白名（Gene Symbol）
organ_pair_ppi_named <- organ_pair_ppi %>%
  left_join(gene_map[, c("X.string_protein_id", "Gene_symbol")], by = c("protein1" = "X.string_protein_id")) %>%
  rename(gene1 = Gene_symbol) %>%
  left_join(gene_map[, c("X.string_protein_id", "Gene_symbol")], by = c("protein2" = "X.string_protein_id")) %>%
  rename(gene2 = Gene_symbol)

# 提取所有器官对的蛋白详情
ppi_cross_organ %>%
  mutate(organ_pair = paste0(pmin(organ_from, organ_to), "_", pmax(organ_from, organ_to))) %>%
  group_by(organ_pair) %>%
  summarise(
    interactions = n(),
    protein_pairs = paste0(protein1, "-", protein2, collapse = "; ")
  ) -> organ_pair_summary

organ_pair_ppi_named <- ppi_cross_organ %>%
  left_join(gene_map[, c("X.string_protein_id", "Gene_symbol")], by = c("protein1" = "X.string_protein_id")) %>%
  rename(gene1 = Gene_symbol) %>%
  left_join(gene_map[, c("X.string_protein_id", "Gene_symbol")], by = c("protein2" = "X.string_protein_id")) %>%
  rename(gene2 = Gene_symbol)



# 器官交互图 ----
# 获取器官名称
ppi_summary <- ppi_cross_organ %>%
  count(organ_from, organ_to, name = "interaction_count") %>%
  mutate(pair = paste0(pmin(organ_from, organ_to), "_", pmax(organ_from, organ_to)))

organ_order <- unique(c(ppi_summary$organ_from, ppi_summary$organ_to))  # 或者用你自己整理好的器官顺序

# 生成足够的 pastel 色
n_organs <- length(organ_order)
pastel_palette <- colorRampPalette(brewer.pal(9, "Pastel1"))(n_organs)

# 分配颜色
organ_colors <- setNames(pastel_palette, organ_order)

# 绘制
circos.clear()

chordDiagram(
  ppi_summary,
  grid.col = organ_colors,
  transparency = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = 1
)
circos.trackPlotRegion(track.index = 1,panel.fun = function(x, y) {
  organ = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + .1, organ, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)

tiff("organ_chord_diagram.tiff", width = 8, height = 8, units = "in", res = 600, compression = "lzw")
chordDiagram(ppi_summary, grid.col = organ_colors, transparency = 0.25,
             annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  organ = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + .1, organ, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)
dev.off()

# Hub器官 ----
top_organs <- organ_centrality %>%
  arrange(desc(Strength)) %>%
  slice(1:10) %>%
  pull(Organ)

subgraph <- induced_subgraph(organ_graph, vids = top_organs)

# 创建 tidygraph 对象
tg <- as_tbl_graph(subgraph)

# 添加中心性指标
tg <- tg %>%
  mutate(Strength = strength(subgraph, mode = "all", weights = E(subgraph)$interaction_count))

# 绘制网络图
ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = interaction_count), color = "gray70", alpha = 0.8) +
  geom_node_point(aes(size = Strength, color = Strength)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, color = "black") +
  scale_size_continuous(range = c(5, 15)) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  theme_void() +
  labs(title = "Hub Organ Network")

# Organ 条形图
# 假设 organ_centrality 已经包含 Organ 和 Strength 两列
# 取 Top10 Hub 器官
df_bar <- organ_centrality %>%
  arrange(desc(Strength)) %>%
  slice(1:10)

ggplot(df_bar, aes(x = reorder(Organ, Strength), y = Strength, fill = Strength)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = round(Strength, 1)), 
            hjust = -0.1, 
            size = 3.5, 
            color = "gray20") +
  scale_fill_gradient(low = "#657DAF", high = "#A7687D") +
  coord_flip() +
  labs(x = NULL, 
       y = "Interaction Strength",
       title = "Top 10 Hub Organs by Network Strength") +
  theme_minimal(base_family = "Arial") +
  theme(
    plot.title     = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text      = element_text(size = 12, color = "gray30"),
    axis.title.y   = element_blank(),
    axis.title.x   = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  expand_limits(y = max(df_bar$Strength) * 1.1)

# 径向条形图
# 仍然取 Top10，并计算角度与对齐方式
df_radial <- df_bar %>%
  mutate(
    Organ = factor(Organ, levels = Organ),
    angle = 90 - 360 * (as.numeric(Organ) - 0.5) / n(),
    hjust  = ifelse(angle < -90, 1, 0)
  )

ggplot(df_radial, aes(x = Organ, y = Strength, fill = Strength)) +
  geom_col(width = 1, show.legend = FALSE) +
  geom_text(aes(label = round(Strength,1), 
                y = Strength + max(Strength)*0.05,
                angle = angle, 
                hjust = hjust),
            size = 3, 
            color = "gray20") +
  scale_fill_gradient(low = "#FDBF6F", high = "#E31A1C") +
  coord_polar(start = 0) +
  labs(x = NULL, y = NULL, title = "Radial Bar Chart of Hub Organs") +
  theme(
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin     = margin(5, 5, 5, 5)
  )




# cell type解析 ----
# 读取Enrichr结果（示例）
enrichr_results <- read.delim("./01_data/Human_Gene_Atlas_table.txt")

# 筛选显著条目（FDR < 0.05）
significant_celltypes <- enrichr_results %>%
  filter(P.value < 0.05) %>%
  arrange(Combined.Score)

# 查看Top10细胞类型
head(significant_celltypes, 10)

library(ggplot2)

ggplot(significant_celltypes[1:12, ], 
       aes(x = Combined.Score, 
           y = reorder(Term, Combined.Score),
           size = -log10(P.value),
           color = Combined.Score)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Enrichment of cell types for aging markers",
       x = "Combined Score", 
       y = "Cell Type") +
  theme_minimal()


## Cell type & Secreted重新解析 ----
# 准备基因列表
# 假设你已经有 overlap_aging_Plasma_secreted 这个向量或数据框
# 转换为数据框，并设定列名

overlap_aging_Plasma_secreted <- overlap_aging_Plasma_secreted %>%
  rename(Genes = 1)
genes <- unique(overlap_aging_Plasma_secreted$Genes)
# cell type注释
celltype_map <- hpa_tissue %>%
  filter(`Gene.name` %in% genes, Level %in% c("High", "Medium")) %>%
  select(Gene = `Gene.name`, Tissue, `Cell.type`, Level) %>%
  distinct()

# secreted注释
secreted_map <- hpa_secretome %>%
  filter(Gene %in% genes & !is.na(Secretome.location)) %>%
  select(Gene, Secretome.location)

# 合并注释信息
gene_annotation <- celltype_map %>%
  left_join(secreted_map, by = "Gene") %>%
  distinct()

# 准备环状图数据框
circle_data <- gene_annotation %>%
  count(`Cell.type`, Secretome.location) %>%
  filter(!is.na(Secretome.location), !is.na(`Cell.type`))


# 按 Cell type 计数
celltype_count <- circle_data %>%
  count(`Cell.type`) %>%
  arrange(desc(n))

# 绘图：横向条形图
ggplot(celltype_count, aes(x = reorder(`Cell.type`, n), y = n)) +
  geom_col(fill = "#4E79A7") +
  coord_flip() +
  labs(title = "All overlap protein Cell Type ",
       x = "Cell type",
       y = "Number") +
  theme_minimal()
ggsave("all overlap protein cell type.pdf", width = 6, height = 12)

### 只保留前10个 Cell type----
celltype_count_top10 <- celltype_count %>%
  slice_max(n, n = 10)

# 绘图
ggplot(celltype_count_top10, aes(x = reorder(`Cell.type`, n), y = n)) +
  geom_col(fill = "#4E79A7") +
  coord_flip() +
  labs(title = "Top 10 Cell Types of Overlap Proteins",
       x = "Cell type",
       y = "Number") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

ggsave("top10_overlap_protein_cell_type.pdf", width = 6, height = 9)


# 分泌蛋白的celltype x secreted locataion关系的桑基图
# 创建冷色调色板（自定义或 RColorBrewer）
n_colors <- length(unique(circle_data$`Cell.type`))
cool_colors <- colorRampPalette(brewer.pal(9, "BuGn"))(n_colors)

# 绘图
ggplot(circle_data,
       aes(axis1 = `Cell.type`, axis2 = Secretome.location)) +
  geom_alluvium(aes(fill = `Cell.type`), width = 0.12, alpha = 0.8) +
  geom_stratum(width = 0.12, fill = "grey95", color = "grey50") +
  geom_text(stat = "stratum",
            aes(label = stringr::str_wrap(after_stat(stratum), width = 18)), 
            size = 3, color = "black", lineheight = 0.9) +
  scale_fill_manual(values = cool_colors) +
  scale_x_discrete(limits = c("Cell type", "Secretome location"),
                   expand = c(0.1, 0.1),
                   name = NULL) +
  labs(title = "Cell Type & Secretome Association",
       y = "Protein number") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
  )
ggsave("Cell Type & Secretome Association.pdf", width = 6, height = 18)

# dot plot
library(ggplot2)
circle_data %>%
  count(`Cell.type`, Secretome.location) %>%
  ggplot(aes(x = Secretome.location, y = `Cell.type`, size = n)) +
  geom_point(alpha = 0.8, color = "#1B9E77") +
  scale_size_continuous(name = "蛋白数量", range = c(2, 10)) +
  theme_classic(base_size = 13) +
  labs(
    title = "不同 Cell Type 分泌到不同位置的蛋白数量分布",
    x = "分泌定位", y = "Cell Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
  )

# secreted location 饼状图
# 汇总数据
secretome_summary <- circle_data %>%
  count(Secretome.location) %>%
  arrange(desc(n)) %>%
  mutate(prop = n / sum(n),
         percent_label = paste0(round(prop * 100, 1), "%"))

# 冷色调调色板
n_colors <- nrow(secretome_summary)
cool_colors <- colorRampPalette(brewer.pal(9, "PuBuGn"))(n_colors)

# 饼图
ggplot(secretome_summary, aes(x = "", y = prop, fill = Secretome.location)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = cool_colors) +
  labs(title = "Secretome Location distribution", fill = "Secretome location") +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "right"
  )
ggsave("Secretome Location distribution.pdf", width = 8, height = 6)

# 生成 Secretome.location 分布表格
secretome_table <- circle_data %>%
  count(Secretome.location) %>%
  arrange(desc(n)) %>%
  mutate(比例 = round(100 * n / sum(n), 1)) %>%
  rename(
    `分泌位置 (Secretome.location)` = Secretome.location,
    `蛋白数量` = n
  )

# 查看表格
print(secretome_table)

# 映射多器官衰老共享marker----

# Step 1: 读取你的基因列表（Gene_symbol + organ）
gene_list_shared <- read_excel("./01_data/shared marker.xlsx")
View(gene_list_shared)

# Step 2: 加载 protein.aliases 文件，手动映射 Gene_symbol -> STRING ID
aliases <- read.delim("./01_data/9606.protein.aliases.v11.5.txt.gz", header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

# 过滤常用命名类型，选取最有代表性的映射
aliases_filtered <- aliases %>% filter(source %in% c("Ensembl_HGNC", "Ensembl", "Gene_Name"))
gene_map_shared <- merge(gene_list_shared, aliases_filtered, by.x="Gene_symbol", by.y="alias")
colnames(gene_map)[colnames(gene_map_shared) == "stringId"] <- "STRING_id"

# Step 3: 加载 protein.links 文件，提取我们需要的 PPI 关系
links <- read.delim("./01_data/9606.protein.links.v11.5.txt.gz", header=TRUE, sep=" ")

# 只保留我们的输入基因之间的交互
input_ids_shared <- gene_map_shared$X.string_protein_id
ppi_filtered_shared <- links %>%
  filter(protein1 %in% input_ids_shared & protein2 %in% input_ids_shared)

# Step 4: 给交互表加上器官信息
ppi_annotated_shared <- ppi_filtered_shared %>%
  left_join(gene_map_shared[, c("X.string_protein_id", "Organ")], by = c("protein1" = "X.string_protein_id")) %>%
  rename(organ_from = Organ) %>%
  left_join(gene_map_shared[, c("X.string_protein_id", "Organ")], by = c("protein2" = "X.string_protein_id")) %>%
  rename(organ_to = Organ)

# Step 5: 提取跨器官的交互
ppi_cross_organ_shared <- ppi_annotated_shared %>% filter(organ_from != organ_to)

# Step 6: 构建器官-器官连接表
organ_edges_shared <- ppi_cross_organ_shared %>%
  group_by(organ_from, organ_to) %>%
  summarise(interaction_count = n()) %>%
  ungroup()

# Step 7: 构建器官网络图
organ_graph_shared <- graph_from_data_frame(organ_edges_shared, directed = FALSE)

# Step 8: 可视化
plot(
  organ_graph_shared,
  vertex.size = 30,
  vertex.label.cex = 1.2,
  edge.width = E(organ_graph_shared)$interaction_count / max(E(organ_graph_shared)$interaction_count) * 10,
  main = "Shared marker Interaction Network via PPI"
)

set.seed(123)
plot(
  organ_graph_shared,
  layout = layout_with_fr(organ_graph_shared),  # Fruchterman-Reingold 布局
  vertex.size = scaled_strength,
  vertex.label.cex = 1.2,
  vertex.label.color = "black",
  vertex.color = "skyblue",
  edge.width = scales::rescale(E(organ_graph_shared)$interaction_count, to = c(1, 10)),
  edge.color = "gray70",
  main = "Shared marker PPI Network (Node size = Strength)"
)

# 中心性分析
# Degree (连接了多少其他器官)
deg_shared <- igraph::degree(organ_graph_shared, mode = "all")

# Strength (总边权和 = 所有交互蛋白数)
strength_shared <- strength(organ_graph_shared, mode = "all", weights = E(organ_graph_shared)$interaction_count)

# Betweenness
btw_shared <- betweenness(organ_graph_shared, weights = 1/E(organ_graph_shared)$interaction_count, normalized = TRUE)

# Closeness
cls_shared <- closeness(organ_graph_shared, weights = 1/E(organ_graph_shared)$interaction_count, normalized = TRUE)

# Eigenvector
eig_shared <- eigen_centrality(organ_graph_shared, weights = E(organ_graph_shared)$interaction_count)$vector

# 组合成表格
organ_centrality_shared <- data.frame(
  Organ = names(deg_shared),
  Degree = deg_shared,
  Strength = strength_shared,
  Betweenness = btw_shared,
  Closeness = cls_shared,
  Eigenvector = eig_shared
)

# 再次绘图
# 添加中心性信息
organ_graph_shared <- set_vertex_attr(organ_graph_shared, "Strength", value = strength_shared)

# 转换为 tidygraph 对象
tg_shared <- as_tbl_graph(organ_graph_shared)
# 创建颜色映射
strength_range_shared <- range(tg_shared %>% activate(nodes) %>% pull(Strength))
tg_shared <- tg_shared %>% 
  mutate(
    strength_scaled_shared = (Strength - strength_range_shared[1]) / diff(strength_range_shared),
    label = name
  )

# 开始绘图
ggraph(tg_shared, layout = "fr") + 
  geom_edge_link(aes(width = interaction_count), color = "gray80", alpha = 0.7) +
  geom_node_point(aes(size = Strength, color = strength_scaled_shared)) +
  geom_node_text(aes(label = label), repel = TRUE, size = 4, color = "black") +
  scale_color_gradientn(colors = brewer.pal(9, "YlOrRd")) +
  scale_size(range = c(5, 15)) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  ggtitle("Organ-level PPI Network of Shared Aging Markers")



#看器官对之间的所有交互蛋白
# 设置器官对
organA <- "Heart"
organB <- "Liver"

# 提取两者之间的交互（注意双向）
organ_pair_ppi_shared <- ppi_cross_organ_shared %>%
  filter((organ_from == organA & organ_to == organB) |
           (organ_from == organB & organ_to == organA))

# 加回蛋白名（Gene Symbol）
organ_pair_ppi_named_shared <- organ_pair_ppi_shared %>%
  left_join(gene_map[, c("X.string_protein_id", "Gene_symbol")], by = c("protein1" = "X.string_protein_id")) %>%
  rename(gene1 = Gene_symbol) %>%
  left_join(gene_map[, c("X.string_protein_id", "Gene_symbol")], by = c("protein2" = "X.string_protein_id")) %>%
  rename(gene2 = Gene_symbol)

# 提取所有器官对的蛋白详情
ppi_cross_organ_shared %>%
  mutate(organ_pair = paste0(pmin(organ_from, organ_to), "_", pmax(organ_from, organ_to))) %>%
  group_by(organ_pair) %>%
  summarise(
    interactions = n(),
    protein_pairs = paste0(protein1, "-", protein2, collapse = "; ")
  ) -> organ_pair_summary_shared

organ_pair_ppi_named_shared <- ppi_cross_organ_shared %>%
  left_join(gene_map[, c("X.string_protein_id", "Gene_symbol")], by = c("protein1" = "X.string_protein_id")) %>%
  rename(gene1 = Gene_symbol) %>%
  left_join(gene_map[, c("X.string_protein_id", "Gene_symbol")], by = c("protein2" = "X.string_protein_id")) %>%
  rename(gene2 = Gene_symbol)


ppi_summary_shared <- ppi_cross_organ_shared %>%
  mutate(organ_pair = paste(pmin(organ_from, organ_to), pmax(organ_from, organ_to), sep = "_")) %>%
  group_by(organ_from, organ_to) %>%
  summarise(interaction_count = n(), .groups = "drop")


# 器官交互图 ----
library(circlize)

# 获取器官名称
organ_order_shared <- unique(c(ppi_summary_shared$organ_from, ppi_summary_shared$organ_to))  # 或者用你自己整理好的器官顺序

# 生成足够的 pastel 色
n_organs_shared <- length(organ_order_shared)
pastel_palette_shared <- colorRampPalette(brewer.pal(9, "Pastel1"))(n_organs_shared)

# 分配颜色
organ_colors_shared <- setNames(pastel_palette_shared, organ_order_shared)

# 绘制
circos.clear()

chordDiagram(
  ppi_summary_shared,
  grid.col = organ_colors_shared,
  transparency = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = 1
)
circos.trackPlotRegion(track.index = 1,panel.fun = function(x, y) {
  organ = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + .1, organ, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)

tiff("organ_chord_diagram.tiff", width = 8, height = 8, units = "in", res = 600, compression = "lzw")
chordDiagram(ppi_summary_shared, grid.col = organ_colors_shared, transparency = 0.25,
             annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  organ = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + .1, organ, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)
dev.off()

# 将shared marker映射至原先的弦图中 ----
organ_order <- unique(c(ppi_summary$organ_from, ppi_summary$organ_to))
n_organs <- length(organ_order)
pastel_palette <- colorRampPalette(brewer.pal(9, "Pastel1"))(n_organs)
organ_colors <- setNames(pastel_palette, organ_order)

ppi_summary$pair <- paste0(pmin(ppi_summary$organ_from, ppi_summary$organ_to), "_",
                           pmax(ppi_summary$organ_from, ppi_summary$organ_to))

ppi_summary_shared$pair <- paste0(pmin(ppi_summary_shared$organ_from, ppi_summary_shared$organ_to), "_",
                                  pmax(ppi_summary_shared$organ_from, ppi_summary_shared$organ_to))

ppi_summary$is_shared <- ppi_summary$pair %in% ppi_summary_shared$pair

ppi_summary$link_color <- ifelse(
  ppi_summary$is_shared,
  organ_colors[ppi_summary$organ_from],
  "gray85"
)

library(circlize)
circos.clear()

tiff("highlighted_shared_chord.tiff", width = 8, height = 8, units = "in", res = 600)

chordDiagram(
  x = ppi_summary[, c("organ_from", "organ_to", "interaction_count")],
  grid.col = organ_colors,
  col = ppi_summary$link_color,
  transparency = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  organ = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + .1, organ, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)

dev.off()


## 重新优化
# 设置 organ colors
organ_order <- unique(c(ppi_summary$organ_from, ppi_summary$organ_to))
n_organs <- length(organ_order)
pastel_palette <- colorRampPalette(brewer.pal(9, "Pastel1"))(n_organs)
organ_colors <- setNames(pastel_palette, organ_order)

# 创建 organ 对的标识符
ppi_summary$pair <- paste0(pmin(ppi_summary$organ_from, ppi_summary$organ_to), "_",
                           pmax(ppi_summary$organ_from, ppi_summary$organ_to))
ppi_summary_shared$pair <- paste0(pmin(ppi_summary_shared$organ_from, ppi_summary_shared$organ_to), "_",
                                  pmax(ppi_summary_shared$organ_from, ppi_summary_shared$organ_to))

# 标记 shared vs non-shared
ppi_summary$is_shared <- ppi_summary$pair %in% ppi_summary_shared$pair

# 设置颜色和透明度
ppi_summary$link_color <- ifelse(
  ppi_summary$is_shared,
  organ_colors[ppi_summary$organ_from],  # 使用 organ color 表示 shared
  "gray85"                               # 非 shared 用灰色
)
ppi_summary$link_transparency <- ifelse(
  ppi_summary$is_shared,
  0.15,  # 不透明，更醒目
  0.6    # 更透明，更淡化
)

# 开始画图
circos.clear()

pdf("highlighted_shared_chord.pdf", width = 8, height = 8)

chordDiagram(
  x = ppi_summary[, c("organ_from", "organ_to", "interaction_count")],
  grid.col = organ_colors,
  col = ppi_summary$link_color,
  transparency = ppi_summary$link_transparency,  # 使用每条边自己的透明度
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  organ = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + .1, organ, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)

dev.off()

# 将各器官蛋白进行归一化 ----
# 假设 organ_edges 是 data.frame(organ_from, organ_to, interaction_count)
# 假设 organ_sizes 是 data.frame(organ, n_proteins)

# 1) 计算每个器官的蛋白数
organ_sizes <- gene_map %>%
  distinct(Organ, Gene_symbol) %>%
  count(Organ, name = "n_proteins")

# 2) 合并到 organ_edges
# 根据两器官蛋白数几何平均进行归一化
edges_norm <- organ_edges %>%
  left_join(organ_sizes,    by = c("organ_from" = "Organ")) %>%
  rename(nA = n_proteins) %>%
  left_join(organ_sizes,    by = c("organ_to"   = "Organ")) %>%
  rename(nB = n_proteins) %>%
  # 3) 归一化
  mutate(interaction_norm = interaction_count / sqrt(nA * nB))

# 根据Min-Max或Z-score进行归一化
edges_norm2 <- organ_edges %>%
  mutate(
    # Min–Max 归一化
    interaction_minmax = (interaction_count - min(interaction_count)) /
      (max(interaction_count) - min(interaction_count)),
    # Z‑score
    interaction_z = (interaction_count - mean(interaction_count)) /
      sd(interaction_count)
  )

# 绘图
library(tidygraph)
library(ggraph)

# 构建 tidygraph 对象
tg_norm <- tbl_graph(
  nodes = data.frame(name = unique(c(edges_norm$organ_from, edges_norm$organ_to))),
  edges = edges_norm %>% select(from = organ_from, to = organ_to, weight = interaction_norm),
  directed = FALSE
)

# 在 nodes 上激活，计算加权度（Strength_norm）
tg_norm <- tg_norm %>%
  activate(nodes) %>%
  mutate(
    Strength_norm = scale(centrality_degree(mode = "all", weights = weight))
  )


# 在 edges 上激活，准备可视化
tg_norm %>%
  activate(edges) %>%
  # 可选：将交互数按权重映射到边宽
  mutate(width_norm = scales::rescale(weight, to = c(0.5, 5))) %>%
  activate(nodes) %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(width = width_norm),
                 color = "gray70", alpha = 0.6) +
  geom_node_point(aes(size = Strength_norm),
                  color = "steelblue") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_width_identity() +
  scale_size_continuous(range = c(3, 12)) +
  theme_void() +
  labs(title = "Normalized Organ-Level PPI Network",
       subtitle = "Edge width = interaction_norm; Node size = scaled Strength")

ggraph(tg_norm, layout = "fr") + 
  # 边：宽度按 weight，颜色固定
  geom_edge_link(aes(width = weight), color = "gray80", alpha = 0.6) +
  # 节点：大小和颜色同时映射到 Strength_norm
  geom_node_point(aes(size = Strength_norm, color = Strength_norm)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, color = "black") +
  # 节点颜色渐变（低值浅蓝，高值深蓝）
  scale_color_gradient(low = "#297270", high = "#e66d50") +
  # 节点大小范围
  scale_size_continuous(range = c(3, 12)) +
  theme_void() +
  labs(title = "Normalized Organ PPI Network", color = "Node centrality", size = "Node centrality")

# KEGG&GO分析----
All_overlap_df <- read_excel("D:/Rproject/organ aging biomarker/result/GO&KEGG/All_overlap.xlsx", col_names = FALSE)

# 提取这一列为基因名向量
gene_list <- All_overlap_df[[1]]  # 只有一列，直接提取即可

# 调用你的富集分析函数
run_enrichment_analysis(
  gene = gene_list,
  OrgDb = "Hs",
  dir = paste0(dir_result, "/Brain")
)

shared_marker_df <- read_excel("D:/Rproject/organ aging biomarker/result/GO&KEGG/Shared_marker.xlsx", col_names = FALSE)

# 提取这一列为基因名向量
gene_list <- shared_marker_df[[1]]  # 只有一列，直接提取即可

# 调用你的富集分析函数
run_enrichment_analysis(
  gene = gene_list,
  OrgDb = "Hs",
  dir = paste0(dir_result, "/shared_marker")
)

# AgingAtlas数据对比分析----
tableExport_1 <- read_csv("~/project/OrganAging/01_data/tableExport.csv")
tableExport_2 <- read_csv("~/project/OrganAging/01_data/tableExport (1).csv")
tableExport_3 <- read_csv("~/project/OrganAging/01_data/tableExport (2).csv")
tableExport_4 <- read_csv("~/project/OrganAging/01_data/tableExport (3).csv")
tableExport_5 <- read_csv("~/project/OrganAging/01_data/tableExport (4).csv")
tableExport_6 <- read_csv("~/project/OrganAging/01_data/tableExport (5).csv")

AgingAtlas <- bind_rows(
  tableExport_1, tableExport_2, tableExport_3,
  tableExport_4, tableExport_5, tableExport_6
)

write.csv(AgingAtlas, file="./03_result/AgingAtlas.csv")

# 与收集的biomarker取交集
# 获取交集
overlap_genes <- intersect(AgingAtlas$Symbol, BIOMARKER$Gene)

# 查看前几项交集
head(overlap_genes)

# 获取两个集合
set1 <- AgingAtlas$Symbol
set2 <- BIOMARKER$Gene

# 计算交集
intersect_len <- length(intersect(set1, set2))

# 创建 Venn 图
venn.plot <- draw.pairwise.venn(
  area1 = length(set1),
  area2 = length(set2),
  cross.area = intersect_len,
  category = c("AgingAtlas", "OrganAging"),
  fill = c("#657DAF", "#A7687D"),
  alpha = 0.6,
  cat.pos = c(0, 0), # 类别标签的位置
  cat.dist = c(0.03, 0.03),
  cat.cex = 1.2,
  cex = 1.5,
  fontface = "bold",
  lty = "blank"
)

# 保存图像到文件
pdf("./03_result/Venn_organ_vs_agingatlas.pdf", width = 6, height = 6)
grid.draw(venn.plot)
dev.off()


gene_list <- AgingAtlas$Symbol
dir_result <- "./"
# 调用你的富集分析函数
run_enrichment_analysis(
  gene = gene_list,
  OrgDb = "Hs",
  dir = paste0(dir_result, "./03_result/AgingAtlas")
)


# 根据Gene_Set画饼状图
# 统计每类的数量
# 假设 AgingAtlas 数据框存在，且包含 Symbol 和 Gene_Set 两列
# 1. 构造饼图数据
gene_set_count <- AgingAtlas %>%
  group_by(Gene_Set) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = round(Count / sum(Count) * 100, 1),
         Label = paste0(Gene_Set, ": ", Count, " (", Percentage, "%)"))

# 颜色
n_colors <- nrow(gene_set_count)
fill_colors <- brewer.pal(min(12, max(3, n_colors)), "Set3")

# 饼图
p_pie <- ggplot(gene_set_count, aes(x = "", y = Count, fill = Gene_Set)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(values = fill_colors) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Gene Set Composition in AgingAtlas") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# 右侧图例注释：构造数据框
legend_df <- gene_set_count %>%
  arrange(desc(Count)) %>%
  mutate(y = rev(seq_along(Gene_Set))*0.5)  # 垂直排序

p_legend <- ggplot(legend_df) +
  geom_point(aes(x = 1, y = y, color = Gene_Set), size = 5) +
  geom_text(aes(x = 1.1, y = y, label = Label), hjust = 0, size = 6) +
  scale_color_manual(values = fill_colors) +
  theme_void() +
  theme(legend.position = "none") +
  xlim(1, 3)  # 留出注释空间

# 组合图形：左饼图 + 右注释
final_plot <- plot_grid(p_pie, p_legend, nrow = 1, rel_widths = c(1, 1.2))

# 展示图形
print(final_plot)
ggsave("./04_figure/AgingAtlas_GeneSet_Pie_WithLegend.pdf", final_plot, width = 10, height = 6)


# SenNet细胞衰老标志物数据库分析 ----
SenNet <- read_excel("01_data/SenNet.xlsx")
View(SenNet)

# 与收集的biomarker取交集
# 获取交集
overlap_genes <- intersect(SenNet$Marker, BIOMARKER$Gene)

# 查看前几项交集
head(overlap_genes)

# 获取两个集合
set1 <- unique(SenNet$Marker)
set2 <- BIOMARKER$Gene

# 计算交集
intersect_len <- length(intersect(set1, set2))

# 创建 Venn 图
venn.plot <- draw.pairwise.venn(
  area1 = length(set1),
  area2 = length(set2),
  cross.area = intersect_len,
  category = c("AgingAtlas", "OrganAging"),
  fill = c("#657DAF", "#A7687D"),
  alpha = 0.6,
  cat.pos = c(0, 0), # 类别标签的位置
  cat.dist = c(0.03, 0.03),
  cat.cex = 1.2,
  cex = 1.5,
  fontface = "bold",
  lty = "blank"
)

# 保存图像到文件
pdf("./04_figure/Venn_organ_vs_SenNet.pdf", width = 6, height = 6)
grid.draw(venn.plot)
dev.off()


gene_list <- SenNet$Marker
dir_result <- "./"
# 调用你的富集分析函数
run_enrichment_analysis(
  gene = gene_list,
  OrgDb = "Hs",
  dir = paste0(dir_result, "./03_result/SenNet")
)

# 统计画图
# 主饼图的数据（按Tissue计数）
# Tissue 频数
tissue_count <- SenNet %>%
  count(Tissue) %>%
  mutate(Percent = n / sum(n) * 100,
         Label = paste0(Tissue, ": ", n, " (", round(Percent, 1), "%)"))

# 低饱和色板
palette1 <- scales::hue_pal(l = 70, c = 40)(nrow(tissue_count))

# 画饼图
ggplot(tissue_count, aes(x = "", y = n, fill = Tissue)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = palette1) +
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "Tissue Composition in SenNet")

# 各Tissue的Sensecence Hallmark的环状饼图
# 统计每个 Tissue 中各个 Hallmark 的数量
# 准备数据（假设你的数据叫 SenNet，列为 Tissue 和 Senescence_Hallmark）
hallmark_count <- SenNet %>%
  count(Tissue, Senescence_Hallmark, name = "Count") %>%
  group_by(Tissue) %>%
  mutate(Percent = Count / sum(Count)) %>%
  ungroup()

# 把 NA 转为因子标签
hallmark_count$Senescence_Hallmark <- fct_na_value_to_level(hallmark_count$Senescence_Hallmark, level = "Unknown")


# 调色板（低饱和度）
palette2 <- c(
  "Cell cycle arrest" = "#cba39b",
  "Cell surface markers" = "#e6c3b3",
  "Changes in morphology" = "#c8d1b6",
  "DNA damage" = "#d8b4a6",
  "DNA damage response" = "#c9dad4",
  "Increased lysosomal content" = "#abc4b3",
  "Metabolic adaptations" = "#dad0aa",
  "Nuclear changes" = "#b9c3d2",
  "Nuclear reorganization" = "#c9ced6",
  "Other" = "#e5d4c0",
  "SASP" = "#b295c5",
  "Upregulation of anti-apoptotic pathways" = "#dec0d4",
  "Unknown" = "#d9d9d9"  # 为 NA 设置颜色
)

# 绘图
hallmark_plot <- ggplot(hallmark_count) +
  geom_arc_bar(
    aes(
      x0 = 0, y0 = 0, r0 = 0.5, r = 1,
      amount = Percent,
      fill = Senescence_Hallmark
    ),
    stat = "pie",
    alpha = 0.9,
    color = "white"
  ) +
  coord_fixed() +
  facet_wrap(~ Tissue, ncol = 4) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines")  # ← 这里加大间隔
  )+
  scale_fill_manual(values = palette2) +
  theme_void(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +
  labs(title = "Senescence Hallmarks in Each Tissue", fill = "Senescence_Hallmark")

print(hallmark_plot)
ggsave("./04_figure/Tissue_Senescence_Hallmark_Donuts.pdf", hallmark_plot, width = 12, height = 8)
