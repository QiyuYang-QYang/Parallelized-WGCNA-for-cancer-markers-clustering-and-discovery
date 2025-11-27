package main

import (
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

const (
	// input1: GCT data path
	gctDataFile = "gene_reads_v10_thyroid.gct.gz"

	// input2: GTF annotation path (length of genes, for TPM)
	gtfAnnotationFile = "gencode.v36.annotation.gtf.gz"

	// output: matrix cleaned
	outputMatrixFile = "clean_thyroid_matrix.csv"

	// Filtering parameter
	// We delete a gene if it's expressed in 90% samples (based on log2(TPM+1) < 1)
	lowExpressionThreshold = 0.9

	// delete genes of low variance percentile
	lowVariancePercentile = 0.25

	//soft threshold: beta
	softPowerBeta = 6.0
)

func main() {
	//PHASE1: Preprocessing the data (parsing & filtering)
	log.Println("Phase 1: Preprocessing the data (parsing & filtering)")
	// parsing GTF annotations (we need gene length for TPM)
	log.Println("Parsing GTF annotation...")
	// We need to perform **streaming parsing** of GTF because it becomes extremely large after decompression.
	// With parseGTFtoLengths, we get: 
	// map[gene_id_with_version] -> length_in_kilobases
	geneLengthsKB, err := parseGTFToLengths(gtfAnnotationFile)
	if err != nil {
		log.Fatalf("Failed: %v", err)
	}
	log.Printf("...succeed in parsing %d genes length\n", len(geneLengthsKB))

	
	// preprocessing GCT raw main counts
	
	log.Println("Preprocessing GCT raw counts")
	
	// With processGCTFile, we get a cleaned matrix.
	finalMatrix, finalGeneList, finalSampleList, err := processGCTFile(
		gctDataFile,
		geneLengthsKB,
		lowExpressionThreshold,
		lowVariancePercentile,
	)
	if err != nil {
		log.Fatalf("Failed: %v", err)
	}

	
	log.Println("Phase 2: Correlation matrix & Adjacency matrix")
	// PHASE2: Pearsons matrix; Correlation matrix
	// ---------------------------------------------------------
	log.Println("Generating the matrix:", outputMatrixFile)
	err = writeOutputCSV(outputMatrixFile, finalMatrix, finalGeneList, finalSampleList)
	if err != nil {
		log.Fatalf("Failed in writing the matrix: %v", err)
	}

	
	log.Printf("  (P2) uses a %d gene x %d sample matrix", len(finalGeneList), len(finalSampleList))
	correlationMatrix, err := RunPhase2(finalMatrix, finalGeneList)
	if err != nil {
		log.Fatalf("failed to run phase 2: %v", err)
	}

	err = writeCorrelationMatrix("correlation_matrix.csv", correlationMatrix, finalGeneList)
	if err != nil {
		// Log a warning if saving fails, but don't stop the program
		log.Printf("warning: failed to save correlation matrix: %v", err)
	}
	// PHASE 3: construct the Adjacency Matrix
	// ---------------------------------------------------------
	log.Println("Phase 3: Calculating Adjacency Matrix...")
    log.Printf(" -> Applying Soft Thresholding with Beta = %.1f", softPowerBeta)

	adjacencyMatrix := CalculateAdjacencyMatrix(correlationMatrix, softPowerBeta)
	log.Printf(" -> Adjacency Matrix created. Size: %d x %d", len(adjacencyMatrix), len(adjacencyMatrix))
	//save the adjacency matrix
	log.Println("Saving Adjacency Matrix to CSV...")
    err = writeCorrelationMatrix("adjacency_matrix.csv", adjacencyMatrix, finalGeneList) 
    if err != nil {
        log.Printf("warning: failed to save adjacency matrix: %v", err)
    }
	// PHASE 4: Topological Overlap Matrix (TOM)
	// ---------------------------------------------------------
	log.Println("Phase 4: Calculating Topological Overlap Matrix (TOM)...")
	tomMatrix := CalculateTOM(adjacencyMatrix)
	log.Printf(" -> TOM created. Size: %d x %d", len(tomMatrix), len(tomMatrix))
	log.Println("Saving TOM Matrix to CSV (This might be large)...")
    err = writeCorrelationMatrix("tom_matrix.csv", tomMatrix, finalGeneList)
    if err != nil {
        log.Printf("warning: failed to save TOM matrix: %v", err)
    }

	// PHASE 5: Prepare for Clustering (Dissimilarity)
    // ---------------------------------------------------------
    log.Println("Phase 5: Calculating Dissimilarity (1-TOM)...")
	distMatrix := CalculateDissimilarity(tomMatrix)
	
	// save Dissimilarity matrix for RShiny visualization.
	finalFile := "dissimilarity_matrix.csv"
	log.Println("Saving Dissimilarity Matrix for clustering...")
	err = writeCorrelationMatrix(finalFile, distMatrix, finalGeneList)
	if err != nil {
		log.Fatalf("Failed to save final dissimilarity matrix: %v", err)
	}

	log.Println("DONE! Pipeline finished.")
}


// writeOutputCSV takes in filepath, processed matrix, genelist, samplelist input
// It saves the matrix in local environment.
func writeOutputCSV(filePath string, matrix [][]float64, geneList []string, sampleList []string) error {
	file, err := os.Create(filePath)
	if err != nil {
		return err
	}
	defer file.Close()

	// Write in sample id
	// add a gene_id column
	header := append([]string{"gene_id"}, sampleList...)
	_, err = fmt.Fprintln(file, strings.Join(header, ","))
	if err != nil {
		return err
	}

	// write in the gene data
	for i, geneID := range geneList {
		row := matrix[i]
		rowStr := make([]string, len(row)+1)
		rowStr[0] = geneID
		for j, val := range row {
			rowStr[j+1] = strconv.FormatFloat(val, 'f', 6, 64) 
		}
		_, err = fmt.Fprintln(file, strings.Join(rowStr, ","))
		if err != nil {
			
			log.Printf("Failed to write %s, %v", geneID, err)
		}
	}
	return nil
}

