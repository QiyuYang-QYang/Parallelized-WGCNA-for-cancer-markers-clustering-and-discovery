package main

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)


const (
	inputRootDir = "./gdc"                             
	outputGCT    = "./gene_reads_v10_thyroid.gct.gz"   
)

func main() {
	log.Println("=== Building GCT from STAR TSV files ===")
	log.Printf("Input directory: %s\n", inputRootDir)
	log.Printf("Output GCT file: %s\n", outputGCT)

	// gene_id → sample → count
	geneCounts := make(map[string]map[string]float64)
	sampleSet := make(map[string]struct{})

	// traverse through all tsv.s in gdc/ 
	err := filepath.Walk(inputRootDir, func(path string, info os.FileInfo, err error) error {
		if err != nil {
			return err
		}
		if info.IsDir() {
			return nil
		}
		if !strings.HasSuffix(info.Name(), ".tsv") {
			return nil
		}

		// sample name
		sampleName := filepath.Base(filepath.Dir(path))
		log.Printf("Parsing sample %s from %s", sampleName, path)
		sampleSet[sampleName] = struct{}{}

		if err := parseOneTSV(path, sampleName, geneCounts); err != nil {
			return fmt.Errorf("failed parsing %s: %w", path, err)
		}
		return nil
	})
	if err != nil {
		log.Fatalf("Walk directory error: %v", err)
	}

	// sample ID and gene ID
	samples := make([]string, 0, len(sampleSet))
	for s := range sampleSet {
		samples = append(samples, s)
	}
	sort.Strings(samples)

	genes := make([]string, 0, len(geneCounts))
	for g := range geneCounts {
		genes = append(genes, g)
	}
	sort.Strings(genes)

	log.Printf("Total samples: %d", len(samples))
	log.Printf("Total genes:   %d", len(genes))

	// write GCT file
	if err := writeGCT(outputGCT, geneCounts, genes, samples); err != nil {
		log.Fatalf("Write GCT failed: %v", err)
	}

	log.Println("🎉 GCT building finished successfully!")
}

// parse single TSV 
func parseOneTSV(path, sample string, geneCounts map[string]map[string]float64) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	r := bufio.NewReader(f)

	// first line may be "# gene-model: ..."
	line, err := r.ReadString('\n')
	if err != nil {
		return err
	}
	if !strings.HasPrefix(line, "#") {

		r = bufio.NewReader(io.MultiReader(strings.NewReader(line), r))
	}

	// separate with tab 
	csvr := csv.NewReader(r)
	csvr.Comma = '\t'
	csvr.ReuseRecord = true

	// header
	header, err := csvr.Read()
	if err != nil {
		return err
	}

	// find gene_id column and unstranded column
	geneIdx := -1
	countIdx := -1
	for i, col := range header {
		if col == "gene_id" {
			geneIdx = i
		}
		if col == "unstranded" {
			countIdx = i
		}
	}

	if geneIdx == -1 || countIdx == -1 {
		return fmt.Errorf("tsv missing gene_id or unstranded column: %s", path)
	}

	// read line by line
	for {
		record, err := csvr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		gid := record[geneIdx]
		if gid == "" {
			continue
		}
		if strings.HasPrefix(gid, "N_") { // N_unmapped / N_multimapping 等跳过
			continue
		}

		valStr := record[countIdx]
		val, err := strconv.ParseFloat(valStr, 64)
		if err != nil {
			continue
		}

		if _, ok := geneCounts[gid]; !ok {
			geneCounts[gid] = make(map[string]float64)
		}
		geneCounts[gid][sample] = val
	}

	return nil
}

// write GCT.gz file
func writeGCT(output string, geneCounts map[string]map[string]float64, genes, samples []string) error {
	f, err := os.Create(output)
	if err != nil {
		return err
	}
	defer f.Close()

	gzw := gzip.NewWriter(f)
	defer gzw.Close()

	w := bufio.NewWriter(gzw)
	defer w.Flush()

	// GCT header
	fmt.Fprintln(w, "#1.2")
	fmt.Fprintf(w, "%d\t%d\n", len(genes), len(samples))

	// column name: Name Description Sample1 Sample2...
	header := []string{"Name", "Description"}
	header = append(header, samples...)
	fmt.Fprintln(w, strings.Join(header, "\t"))

	// dataa
	for _, gid := range genes {
		row := []string{gid, gid}
		for _, s := range samples {
			v := geneCounts[gid][s]
			row = append(row, strconv.FormatFloat(v, 'f', 0, 64))
		}
		fmt.Fprintln(w, strings.Join(row, "\t"))
	}

	return nil
}
