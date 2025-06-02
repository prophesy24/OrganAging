library(readxl)
BIOMARKER <- read_excel("rawdata/BIOMARKER.xlsx")
View(BIOMARKER)
SXF_protein <- read_excel("rawdata/SXF protein.xlsx")
View(SXF_protein)
BIOMARKER <- as.data.frame(BIOMARKER)
SXF_protein <- as.data.frame(SXF_protein)

# 提取基因列表并去重
organ_aging <- unique(BIOMARKER$Genes)
SXF_protein <- unique(SXF_protein$Genes)

# 获取共有基因
common_genes <- intersect(organ_aging, SXF_protein)
print(common_genes)

common_genes <- as.data.frame(common_genes)
# 安装并加载包
install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)
# 生成Venn图
ggVennDiagram(
  list(Dataset1 = organ_aging, Dataset2 = SXF_protein),
  category.names = c("Organ aging biomarker", "Plasma protein"),
  set_color = "black",
  set_size = 5,
  label_color = "firebrick",
  edge_size = 1
) +
  labs(title = "Gene Overlap Analysis") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.background = element_rect(fill = "white")
  ) +
  scale_fill_distiller(palette = "RdYlBu")

ggsave("ggVenn.png", width = 6, height = 6, dpi = 300)

# 输出共有基因
write.csv(data.frame(Common_Genes = common_genes), 
          "common_genes.csv", 
          row.names = FALSE)

# 与HPA的分泌蛋白overlap ----

library(hpar)
# 获取HPA所有蛋白质分类数据（含分泌蛋白标记）
secretome <- read.table("./rawdata/sa_location_Secreted.tsv", header=TRUE, sep="\t")

# 加载包
library(dplyr)



# 取交集：筛选你的标志物中属于分泌蛋白的基因
overlap_genes <- inner_join(BIOMARKER, secretome, by = "Gene")

# 保存结果
write.csv(overlap_genes, "secretome_biomarkers_overlap.csv", row.names = FALSE)

# 统计结果
cat(sprintf("你的标志物中有%d个分泌蛋白（占比%.1f%%）",
            nrow(overlap_genes),
            nrow(overlap_genes)/nrow(BIOMARKER)*100))

library(ggplot2)

# 绘制分泌蛋白占比饼图
# 优化后的饼图
ggplot(df, aes(x = "", y = Count, fill = Category)) +
  +     geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +  # 添加白色边框
  +     coord_polar("y", start = 0) +
  +     geom_text(aes(label = Label), 
                  +               position = position_stack(vjust = 0.5),  # 标签居中
                  +               color = "white",  # 标签文字颜色
                  +               size = 5,         # 标签字号
                  +               fontface = "bold") + 
  +     scale_fill_manual(values = my_colors) +  # 应用自定义颜色
  +     labs(title = "Secretome Biomarkers in Organ aging Signature",
             +          subtitle = paste("Total biomarkers:", nrow(BIOMARKER)),
             +          fill = "Protein Category") +
  +     theme_void() +
  +     theme(
    +         plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中加粗
    +         plot.subtitle = element_text(hjust = 0.5, size = 10),  # 副标题
    +         legend.position = "bottom",  # 图例位置
    +         legend.title = element_text(face = "bold")  # 图例标题加粗
    +     )

common_genes$group <- "plasma aging marker"
overlap_genes$group <- "secreted aging marker"

# combined
venn_pre <- rbind(common_genes,overlap_genes)

# res output
write.xlsx(venn_pre, file = "./result/plasma&secreted.xlsx")

secreted <- secretome$Gene
# 分泌蛋白和marker进行overlap
organ_aging <- as.data.frame(organ_aging)
secreted <- as.data.frame(secreted)
colnames(secreted) <- c("Gene")
colnames(organ_aging) <- c("Gene")
organ_aging$group <- "Organ aging marker"
secreted$group <- "Secreted marker"

# combined
venn_pre <- rbind(organ_aging,secreted)
overlap_aging_secreted <- intersect(organ_aging$Gene, secreted$Gene)
length(overlap_aging_secreted)
# res output
write.xlsx(venn_pre, file = "./result/organ aging&secreted.xlsx")

# 血浆蛋白和marker进行overlap
SXF_protein <- as.data.frame(SXF_protein)
colnames(SXF_protein) <- c("Gene")
SXF_protein$group <- "Plasma protein"
# combined
venn_pre <- rbind(organ_aging,SXF_protein)
overlap_aging_Plasma <- intersect(organ_aging$Gene, SXF_protein$Gene)
length(overlap_aging_Plasma)

# res output
write.xlsx(venn_pre, file = "./result/organ aging&plasma.xlsx")

# 将全部overlap
# combined
venn_pre <- rbind(organ_aging,SXF_protein,secreted)
overlap_aging_Plasma_secreted <- intersect(overlap_aging_Plasma, overlap_aging_secreted)
# res output
overlap_aging_Plasma_secreted <- data.frame(overlap_aging_Plasma_secreted)
colnames(overlap_aging_Plasma_secreted) <- overlap_aging_Plasma_secreted$id
write.xlsx(overlap_aging_Plasma_secreted, file = "./result/All overlap.xlsx")

#将基因映射到器官 ----
# 安装并加载hpar包
if (!require("hpar")) BiocManager::install("hpar")
library(hpar)
#allHparData()
#Title
#1          hpaCancer16.1
#2              hpaCancer
#3    hpaNormalTissue16.1
#4        hpaNormalTissue
#5           hpaSecretome
#6    hpaSubcellularLoc14
#7  hpaSubcellularLoc16.1
#8      hpaSubcellularLoc
#9     rnaConsensusTissue
#10       rnaFantomTissue
#11   rnaGeneCellLine16.1
#12       rnaGeneCellLine
#13     rnaGeneTissue21.0
#14         rnaGtexTissue
#15          rnaHpaTissue
# 获取HPA正常组织表达数据（含分泌蛋白信息）
# 使用getHpa()函数获取分泌蛋白数据
#组织数据
download.file("https://www.proteinatlas.org/download/normal_tissue.tsv.zip", "tissue.zip")
unzip("tissue.zip")
hpa_tissue <- read.delim("normal_tissue.tsv")

# 分泌蛋白数据
download.file("https://www.proteinatlas.org/download/proteinatlas.tsv.zip", "secretome.zip")
unzip("secretome.zip")
hpa_secretome <- read.delim("proteinatlas.tsv") %>%
  filter(Secretory == "yes")
# 步骤2：读取数据
library(readr)
hpa_tissue <- read_tsv("normal_ihc_data.tsv")
hpa_secretome <- secretome
# 验证列名
colnames(hpa_tissue)      # 应包含 Gene, Tissue, Cell.type, Level
colnames(hpa_secretome)   # 应包含 Gene, Secretory.score

library(ggplot2)
library(tibble)  # 提供column_to_rownames()
library(dplyr)   # 提供数据操作

# 统计不同分泌定位的蛋白数量
location_count <- hpa_secretome %>%
  filter(!is.na(Secretome.location)) %>%
  count(Secretome.location) %>%
  arrange(desc(n))

# 绘制柱状图
ggplot(location_count, aes(x = reorder(Secretome.location, n), y = n)) +
  geom_col(fill = "#4E79A7") +
  labs(title = "分泌蛋白的亚定位分布", 
       x = "分泌定位", 
       y = "蛋白数量") +
  coord_flip() +  # 横向显示
  theme_minimal()

library(pheatmap)

# 准备数据：筛选高表达分泌蛋白的组织分布
heatmap_data <- hpa_tissue %>%
  filter(`Gene name` %in% secreted_genes & Level %in% c("High", "Medium")) %>%
  count(`Gene name`, Tissue) %>%
  tidyr::pivot_wider(names_from = Tissue, values_from = n, values_fill = 0) %>%
  column_to_rownames("Gene name")

# 绘制热图
pheatmap(heatmap_data, 
         color = colorRampPalette(c("white", "#E15759"))(100),
         main = "分泌蛋白的组织分布热图",
         cluster_rows = TRUE,
         cluster_cols = TRUE)

library(ggraph)
library(tidygraph)

# 构建边数据
edges <- hpa_tissue %>%
  filter(`Gene name` %in% secreted_genes & Level %in% c("High", "Medium")) %>%
  select(from = `Gene name`, to = Tissue)

# 创建网络图
graph <- as_tbl_graph(edges) %>%
  mutate(degree = centrality_degree())

# 绘制网络
set.seed(123)
ggraph(graph, layout = "fr") + 
  geom_edge_link(alpha = 0.2) +
  geom_node_point(aes(size = degree, color = ifelse(name %in% secreted_genes, "Gene", "Tissue"))) +
  geom_node_text(aes(label = ifelse(degree > 10, name, "")), repel = TRUE) +
  scale_color_manual(values = c("#4E79A7", "#F28E2B")) +
  labs(title = "分泌蛋白-器官关联网络") +
  theme_graph()

# 将我的数据映射到器官 ----
# 加载必要的包
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)  # 桑基图/冲击图绘制
library(networkD3)   # 交互式桑基图

# 从HPA加载组织表达数据（若已有hpa_tissue数据框则跳过）
data(hpaNormalTissue)  # 使用hpar包内置数据
# 或从文件读取：
# hpa_tissue <- read.delim("normal_tissue.tsv")

# 筛选高/中表达水平的记录
hpa_filtered <- hpa_tissue %>%
  filter(Level %in% c("High", "Medium"))
hpa_filtered$ENS <- hpa_filtered$Gene
hpa_filtered <- subset(hpa_filtered,select = -c(Gene))
hpa_filtered$Gene <- hpa_filtered$`Gene name`

# 映射基因到器官
gene_organ_mapping <- overlap_genes %>%
  left_join(hpa_filtered, by = "Gene") %>%
  filter(!is.na(Tissue)) %>%  # 移除无组织信息的基因
  distinct(Gene, group, Tissue)  # 去重

# 计算基因-器官-分组的频数
sankey_data <- gene_organ_mapping %>%
  count(group, Tissue, name = "Count")

# 绘制静态桑基图
ggplot(sankey_data,
       aes(axis1 = group,   # 第一层：分组
           axis2 = Tissue,  # 第二层：器官
           y = Count)) +
  geom_alluvium(aes(fill = group), width = 1/12) +  # 流动带
  geom_stratum(width = 1/12, fill = "grey80", color = "grey") +  # 分层块
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +  # 标签
  scale_x_discrete(limits = c("Group", "Tissue"), expand = c(0.05, 0.05)) +
  labs(title = "分泌蛋白分组-器官分布") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 方法2准备节点列表
nodes <- data.frame(
  name = c(unique(sankey_data$group), unique(sankey_data$Tissue)),
  group = c(rep("Group", n_distinct(sankey_data$group)), 
            rep("Tissue", n_distinct(sankey_data$Tissue)))
)

# 准备连接数据
links <- sankey_data %>%
  mutate(
    source = match(group, nodes$name) - 1,  # 从0开始索引
    target = match(Tissue, nodes$name) - 1,
    value = Count
  )

#  绘制交互式桑基图
sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  NodeGroup = "group",
  units = "genes",
  fontSize = 12,
  nodeWidth = 30
)

# 跨器官衰老网络分析 ----
# 安装和加载必要包

library(STRINGdb)
library(igraph)
library(tidyverse)

library(readxl)
library(igraph)
library(tidyverse)

# Step 1: 读取你的基因列表（Gene_symbol + organ）
gene_list <- read_excel("rawdata/The cross-organ aging network.xlsx")
View(gene_list)

# Step 2: 加载 protein.aliases 文件，手动映射 Gene_symbol -> STRING ID
aliases <- read.delim("./rawdata/9606.protein.aliases.v11.5.txt.gz", header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

# 过滤常用命名类型，选取最有代表性的映射
aliases_filtered <- aliases %>% filter(source %in% c("Ensembl_HGNC", "Ensembl", "Gene_Name"))
gene_map <- merge(gene_list, aliases_filtered, by.x="Gene_symbol", by.y="alias")
colnames(gene_map)[colnames(gene_map) == "stringId"] <- "STRING_id"

# Step 3: 加载 protein.links 文件，提取我们需要的 PPI 关系
links <- read.delim("./rawdata/9606.protein.links.v11.5.txt.gz", header=TRUE, sep=" ")

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
deg <- degree(organ_graph, mode = "all")

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
library(tidygraph)
library(ggraph)
library(ggplot2)
library(RColorBrewer)
library(tidygraph)

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

library(dplyr)
library(ggplot2)
library(scales)  # 用于rescale函数

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
library(gganatogram)
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

library(shinybody)

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

library(shinybody)

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
library(htmlwidgets)

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
  library(shinybody)

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
library(circlize)

# 获取器官名称
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
library(ggplot2)

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
library(ggplot2)
library(dplyr)

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
enrichr_results <- read.delim("./rawdata/Human_Gene_Atlas_table.txt")

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
library(dplyr)
# cell type注释
celltype_map <- hpa_tissue %>%
  filter(`Gene name` %in% genes, Level %in% c("High", "Medium")) %>%
  select(Gene = `Gene name`, Tissue, `Cell type`, Level) %>%
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
  count(`Cell type`, Secretome.location) %>%
  filter(!is.na(Secretome.location), !is.na(`Cell type`))


# 按 Cell type 计数
celltype_count <- circle_data %>%
  count(`Cell type`) %>%
  arrange(desc(n))

# 绘图：横向条形图
ggplot(celltype_count, aes(x = reorder(`Cell type`, n), y = n)) +
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
ggplot(celltype_count_top10, aes(x = reorder(`Cell type`, n), y = n)) +
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
library(ggalluvial)
# 创建冷色调色板（自定义或 RColorBrewer）
n_colors <- length(unique(circle_data$`Cell type`))
cool_colors <- colorRampPalette(brewer.pal(9, "BuGn"))(n_colors)

# 绘图
ggplot(circle_data,
       aes(axis1 = `Cell type`, axis2 = Secretome.location)) +
  geom_alluvium(aes(fill = `Cell type`), width = 0.12, alpha = 0.8) +
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
  count(`Cell type`, Secretome.location) %>%
  ggplot(aes(x = Secretome.location, y = `Cell type`, size = n)) +
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
library(dplyr)

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
library(STRINGdb)
library(igraph)
library(tidyverse)

library(readxl)
library(igraph)
library(tidyverse)

# Step 1: 读取你的基因列表（Gene_symbol + organ）
gene_list_shared <- read_excel("rawdata/shared marker.xlsx")
View(gene_list_shared)

# Step 2: 加载 protein.aliases 文件，手动映射 Gene_symbol -> STRING ID
aliases <- read.delim("./rawdata/9606.protein.aliases.v11.5.txt.gz", header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

# 过滤常用命名类型，选取最有代表性的映射
aliases_filtered <- aliases %>% filter(source %in% c("Ensembl_HGNC", "Ensembl", "Gene_Name"))
gene_map_shared <- merge(gene_list_shared, aliases_filtered, by.x="Gene_symbol", by.y="alias")
colnames(gene_map)[colnames(gene_map_shared) == "stringId"] <- "STRING_id"

# Step 3: 加载 protein.links 文件，提取我们需要的 PPI 关系
links <- read.delim("./rawdata/9606.protein.links.v11.5.txt.gz", header=TRUE, sep=" ")

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
deg_shared <- degree(organ_graph_shared, mode = "all")

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
library(tidygraph)
library(ggraph)
library(ggplot2)
library(RColorBrewer)
library(tidygraph)

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
  geom_node_point(aes(size = Strength, color = strength_scaled)) +
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
