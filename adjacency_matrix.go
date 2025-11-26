package main

import (
	"math"
)

func CalculateAdjacencyMatrix(corrMatrix [][]float64, beta float64) [][]float64 {
	numGenes := len(corrMatrix)
	
	
	adjMatrix := make([][]float64, numGenes)
	for i := range adjMatrix {
		adjMatrix[i] = make([]float64, numGenes)
	}

	
	for i := 0; i < numGenes; i++ {
		for j := i; j < numGenes; j++ {
			correlation := corrMatrix[i][j]

			// WGCNA
			weight := math.Pow(math.Abs(correlation), beta)

			adjMatrix[i][j] = weight
			//diagonal
			if i != j {
				adjMatrix[j][i] = weight
			}
		}
	}

	return adjMatrix
}