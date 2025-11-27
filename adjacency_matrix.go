package main

import (
	"math"
)

// CalculateAdjacencyMatrix constructs a *signed* adjacency matrix
// from a Pearson correlation matrix using soft-thresholding power beta.
//
// Signed network：先把相关系数从 [-1,1] 映射到 [0,1]:
//   s_ij = (corr_ij + 1) / 2
// 再做 a_ij = s_ij^beta
//
// 对角线设为 0（WGCNA 里 adjacency 通常 diag=0，后面 TOM 会把 diag 设成 1）。
func CalculateAdjacencyMatrix(corrMatrix [][]float64, beta float64) [][]float64 {
	numGenes := len(corrMatrix)

	adjMatrix := make([][]float64, numGenes)
	for i := range adjMatrix {
		adjMatrix[i] = make([]float64, numGenes)
	}

	for i := 0; i < numGenes; i++ {
		for j := i; j < numGenes; j++ {
			if i == j {
				// 自连边在 adjacency 里设 0，TOM 里再处理成 1
				adjMatrix[i][j] = 0.0
				continue
			}

			corr := corrMatrix[i][j] // ∈ [-1, 1]

			// 映射到 signed 权重：[-1,1] → [0,1]
			s := (corr + 1.0) / 2.0
			if s < 0 {
				s = 0
			} else if s > 1 {
				s = 1
			}

			// soft thresholding
			weight := math.Pow(s, beta)

			adjMatrix[i][j] = weight
			adjMatrix[j][i] = weight
		}
	}

	return adjMatrix
}
