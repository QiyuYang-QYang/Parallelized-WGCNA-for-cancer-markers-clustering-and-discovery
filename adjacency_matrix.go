package main

import (
	"math"
)

// CalculateAdjacencyMatrix constructs a *signed* adjacency matrix
// from a Pearson correlation matrix using soft-thresholding power beta.
//
// Signed network:
// First, map the correlation coefficients from [-1, 1] to [0, 1]:
//   s_ij = (corr_ij + 1) / 2
// Then apply the power operation: a_ij = s_ij^beta
//
// Set the diagonal to 0 (in WGCNA, the adjacency matrix typically has diag=0,
// while the TOM calculation will later set diag to 1).
func CalculateAdjacencyMatrix(corrMatrix [][]float64, beta float64) [][]float64 {
	numGenes := len(corrMatrix)

	adjMatrix := make([][]float64, numGenes)
	for i := range adjMatrix {
		adjMatrix[i] = make([]float64, numGenes)
	}

	for i := 0; i < numGenes; i++ {
		for j := i; j < numGenes; j++ {
			if i == j {
				// later will further adjust in TOM analysis
				adjMatrix[i][j] = 0.0
				continue
			}

			corr := corrMatrix[i][j] // ∈ [-1, 1]

			// signed weight：[-1,1] → [0,1]
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
