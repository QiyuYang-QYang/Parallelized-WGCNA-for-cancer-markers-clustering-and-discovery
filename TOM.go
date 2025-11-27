package main

import (
	"log"
	"math"
	"runtime"
	"sync"
	"time"
)

// CalculateTOM computes the Topological Overlap Matrix (TOM)
// from a symmetric adjacency matrix
// We can construct signed or unsigned, here we choose *signed*
//
//  WGCNA weighted TOM 
//
//   k_i        = sum_{u != i} a_{iu}
//   numerator  = sum_{u != i,j} a_{iu} * a_{ju} + a_{ij}
//   denominator= min(k_i, k_j) + 1 - a_{ij}
//   TOM_{ij}   = numerator / denominator
//
// TOM_{ii} = 1.0
func CalculateTOM(adjMatrix [][]float64) [][]float64 {
	numGenes := len(adjMatrix)
	log.Printf("Starting TOM calculation for %d genes...", numGenes)
	startTime := time.Now()

	// 1. connectivity：k_i = sum_{u != i} a_{iu}
	k := make([]float64, numGenes)
	for i := 0; i < numGenes; i++ {
		sum := 0.0
		for j := 0; j < numGenes; j++ {
			if j == i {
				continue
			}
			sum += adjMatrix[i][j]
		}
		k[i] = sum
	}
	log.Println("...Connectivity (k) calculated.")

	// 2. TOM matrix
	tomMatrix := make([][]float64, numGenes)
	for i := range tomMatrix {
		tomMatrix[i] = make([]float64, numGenes)
	}

	// 3. parallel calculation
	numCPU := runtime.NumCPU()
	log.Printf("...Parallelizing TOM calculation using %d CPUs", numCPU)

	var wg sync.WaitGroup
	progressChan := make(chan int, numGenes)

	worker := func(startRow, endRow int) {
		defer wg.Done()
		for i := startRow; i < endRow; i++ {

			for j := i; j < numGenes; j++ {

				//TOM(i,i) = 1
				if i == j {
					tomMatrix[i][j] = 1.0
					continue
				}

				// numerator = sum_{u != i,j} a_{iu} a_{ju} + a_{ij}
				dotProduct := 0.0
				for u := 0; u < numGenes; u++ {
					if u == i || u == j {
						continue
					}
					dotProduct += adjMatrix[i][u] * adjMatrix[j][u]
				}
				numerator := dotProduct + adjMatrix[i][j]

				// denominator: min(k_i, k_j) + 1 - a_{ij}
				minK := math.Min(k[i], k[j])
				denominator := minK + 1.0 - adjMatrix[i][j]

				// TOM value
				tomValue := 0.0
				if denominator > 0 {
					tomValue = numerator / denominator
				}

				// value =  [0,1]
				if tomValue < 0 {
					tomValue = 0
				} else if tomValue > 1 {
					tomValue = 1
				}

				tomMatrix[i][j] = tomValue
				if i != j {
					tomMatrix[j][i] = tomValue
				}
			}

			progressChan <- 1
		}
	}

	rowsPerWorker := (numGenes + numCPU - 1) / numCPU
	for w := 0; w < numCPU; w++ {
		start := w * rowsPerWorker
		end := start + rowsPerWorker
		if start >= numGenes {
			break
		}
		if end > numGenes {
			end = numGenes
		}
		wg.Add(1)
		go worker(start, end)
	}

	go func() {
		processed := 0
		lastLog := 0
		for range progressChan {
			processed++
			percent := (processed * 100) / numGenes
			if percent >= lastLog+5 {
				log.Printf("...TOM Progress: %d%% (%d/%d rows)", percent, processed, numGenes)
				lastLog = percent
			}
		}
	}()

	wg.Wait()
	close(progressChan)

	duration := time.Since(startTime)
	log.Printf("TOM calculation finished in %v", duration)

	return tomMatrix
}

// CalculateDissimilarity converts TOM into a dissimilarity matrix 
// via dist = 1 - TOM.
func CalculateDissimilarity(tomMatrix [][]float64) [][]float64 {
	rows := len(tomMatrix)
	distMatrix := make([][]float64, rows)
	for i := range distMatrix {
		distMatrix[i] = make([]float64, rows)
		for j := 0; j < rows; j++ {
			distMatrix[i][j] = 1.0 - tomMatrix[i][j]
		}
	}
	return distMatrix
}
