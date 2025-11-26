import pandas as pd
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# 1. 加载 Go 生成的差异矩阵
print("Loading dissimilarity matrix...")
try:
    # index_col=0 意味着第一列是基因名
    df = pd.read_csv("dissimilarity_matrix.csv", index_col=0)
except FileNotFoundError:
    print("Error: dissimilarity_matrix.csv not found. Run the Go program first.")
    sys.exit(1)

gene_names = df.index.tolist()
print(f"Loaded matrix with {len(gene_names)} genes.")

# 2. 层次聚类 (Hierarchical Clustering)
# method='average' 是 WGCNA 的标准做法 (Average Linkage)
# metric='precomputed' 告诉 scipy 我们给的已经是距离矩阵了，不需要它再算欧式距离
print("Running hierarchical clustering...")
linkage_matrix = sch.linkage(df.values, method='average') 
# 注意：如果矩阵非常大，这里可能会报错 "Distance matrix must be symmetric"。
# 如果 Go 计算精度导致微小误差，可以使用: sch.linkage(sch.distance.squareform(df.values), method='average')

# 3. 绘制树状图 (Dendrogram)
print("Plotting dendrogram...")
plt.figure(figsize=(15, 8))
dendrogram = sch.dendrogram(linkage_matrix, labels=gene_names, no_plot=True)

# 我们可以简单画一个概览图
plt.title('Gene Co-expression Clustering Dendrogram')
plt.xlabel('Genes')
plt.ylabel('Dissimilarity (1-TOM)')
sch.dendrogram(linkage_matrix, no_labels=True) # 基因太多时不显示标签
plt.axhline(y=0.9, c='r', ls='--') # 假设的剪切线，用来辅助观察
plt.savefig("dendrogram.png", dpi=300)
print("Dendrogram saved to dendrogram.png")

# 4. 识别模块 (Module Identification)
# 这一步是“剪枝”。我们可以设定一个高度阈值 (t)，或者设定想要的簇数量。
# 0.9 是一个示例阈值，意味着差异度 > 0.9 的才会被分开
cut_height = 0.95 
labels = sch.fcluster(linkage_matrix, t=cut_height, criterion='distance')

# 5. 保存结果
results = pd.DataFrame({
    'GeneID': gene_names,
    'Module_Label': labels
})

# 按模块排序，方便查看
results = results.sort_values(by='Module_Label')
results.to_csv("gene_modules.csv", index=False)

print(f"Clustering done! Identified {results['Module_Label'].nunique()} modules.")
print("Results saved to gene_modules.csv")