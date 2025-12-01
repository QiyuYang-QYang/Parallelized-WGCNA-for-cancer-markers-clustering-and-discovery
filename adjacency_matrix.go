package main

import (
	"math"
)

// CalculateAdjacencyMatrix constructs a *signed* adjacency matrix
// from a Pearson correlation matrix using soft-thresholding power beta.
//
// CORRECTED VERSION - Based on WGCNA official paper:
// Langfelder & Horvath (2008) BMC Bioinformatics 9:559
//
// Signed network formula:
//   s_ij = (1 + cor_ij) / 2    [transforms [-1,1] to [0,1]]
//   a_ij = s_ij^beta
//
// CRITICAL FIX: Diagonal elements must be 1.0 (not 0.0)
func CalculateAdjacencyMatrix(corrMatrix [][]float64, beta float64) [][]float64 {

    numGenes := len(corrMatrix)

    adjMatrix := make([][]float64, numGenes)
    for i := range adjMatrix {
        adjMatrix[i] = make([]float64, numGenes)
    }

    for i := 0; i < numGenes; i++ {
        for j := i; j < numGenes; j++ {

            if i == j {
                adjMatrix[i][j] = 1.0
                continue
            }

            corr := corrMatrix[i][j]

            var weight float64
            if corr > 0 {
                weight = math.Pow(corr, beta)
            } else {
                weight = 0.0
            }

            adjMatrix[i][j] = weight
            adjMatrix[j][i] = weight
        }
    }
    return adjMatrix
}

// CalculateAdjacencyMatrix_Unsigned creates unsigned network
// Formula: a_ij = |cor_ij|^beta
func CalculateAdjacencyMatrix_Unsigned(corrMatrix [][]float64, beta float64) [][]float64 {
	numGenes := len(corrMatrix)

	adjMatrix := make([][]float64, numGenes)
	for i := range adjMatrix {
		adjMatrix[i] = make([]float64, numGenes)
	}

	for i := 0; i < numGenes; i++ {
		for j := i; j < numGenes; j++ {
			if i == j {
				// Diagonal is 1
				adjMatrix[i][j] = 1.0
				continue
			}

			// Unsigned: take absolute value
			absCorr := math.Abs(corrMatrix[i][j])
			
			// Apply power
			weight := math.Pow(absCorr, beta)

			adjMatrix[i][j] = weight
			adjMatrix[j][i] = weight
		}
	}

	return adjMatrix
}