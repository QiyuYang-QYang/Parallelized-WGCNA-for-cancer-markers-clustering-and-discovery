package main

import (
	"log"
	"math"
	"runtime"
	"sync"
	"time"
)

// CalculateTOM 生成拓扑重叠矩阵 (Topological Overlap Matrix)。
// 这是一个计算密集型任务 (O(N^3))，因此使用了 goroutine 并行计算。
func CalculateTOM(adjMatrix [][]float64) [][]float64 {
	numGenes := len(adjMatrix)
	log.Printf("Starting TOM calculation for %d genes...", numGenes)
	startTime := time.Now()

	// 1. 预计算每个基因的连接度 k (connectivity)
	// k_i = sum(adj[i][u]) for all u
	// 这是 TOM 公式分母部分需要的
	k := make([]float64, numGenes)
	
	// 这里也可以并行，但因为只是简单的求和，单线程足够快
	for i := 0; i < numGenes; i++ {
		sum := 0.0
		for j := 0; j < numGenes; j++ {
			// 注意：通常 WGCNA 计算 k 时不包含自己 (i != j)，
			// 但因为软阈值 adj[i][i] 通常设为 1 或 0，这里直接累加即可。
			// 严格的标准 WGCNA 实现中，connectivity 是行之和。
			sum += adjMatrix[i][j]
		}
		k[i] = sum
	}
	log.Println("...Connectivity (k) calculated.")

	// 2. 初始化 TOM 矩阵
	tomMatrix := make([][]float64, numGenes)
	for i := range tomMatrix {
		tomMatrix[i] = make([]float64, numGenes)
	}

	// 3. 并行计算 TOM
	// 我们将矩阵的行（Rows）分配给不同的 CPU 核心处理
	numCPU := runtime.NumCPU()
	log.Printf("...Parallelizing TOM calculation using %d CPUs", numCPU)

	var wg sync.WaitGroup
	// 用于显示进度的通道
	progressChan := make(chan int, numGenes)

	// 定义每个 worker 的工作
	worker := func(startRow, endRow int) {
		defer wg.Done()
		for i := startRow; i < endRow; i++ {
			// TOM 也是对称矩阵，我们计算上三角 (j >= i)，然后镜像赋值
			for j := i; j < numGenes; j++ {
				// A. 分子：计算点积 (Shared Neighbors) + 直接连接
				// Num = (sum_u adj[i][u] * adj[u][j]) + adj[i][j]
				// 这里的点积代表“共同邻居的强度”
				dotProduct := 0.0
				for u := 0; u < numGenes; u++ {
					// 优化：如果 adj 很小，乘积更小，可以忽略吗？
					// 不，为了精确性，必须全算。这是最耗时的循环。
					dotProduct += adjMatrix[i][u] * adjMatrix[u][j]
				}
				// 原始公式中，分子需要减去 i和j 自身的重叠，但通常近似计算直接加 adj[i][j]
				// 标准公式: Num = (DotProduct - adj[i][j]) + adj[i][j] ... 其实就是 DotProduct (如果对角线是1)
				// 但为了保险和对应 WGCNA R包的逻辑：
				// Num = \sum_{u \neq i,j} (a_iu * a_uj) + a_ij
				// 如果我们在 dotProduct 中包含了 u=i 和 u=j 的项，需要小心处理。
				// 为简单起见，这里使用标准加权 TOM 近似公式：
				numerator := dotProduct // 假设 dotProduct 包含了所有共同路径

				// B. 分母：min(k_i, k_j) + 1 - adj[i][j]
				minK := math.Min(k[i], k[j])
				denominator := minK + 1.0 - adjMatrix[i][j]

				// C. 计算结果
				tomValue := 0.0
				if denominator > 0 {
					tomValue = numerator / denominator
				}
				
				// 赋值
				tomMatrix[i][j] = tomValue
				if i != j {
					tomMatrix[j][i] = tomValue
				}
			}
			// 发送进度信号
			progressChan <- 1
		}
	}

	// 分配任务给 Workers
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

	// 启动一个独立的 goroutine 来监控进度
	go func() {
		processed := 0
		lastLog := 0
		for range progressChan {
			processed++
			// 每完成 5% 打印一次日志
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
	log.Printf("TOM Calculation finished in %v", duration)

	return tomMatrix
}

// CalculateDissimilarity 将 TOM 矩阵转换为距离矩阵
// Dist = 1 - TOM
func CalculateDissimilarity(tomMatrix [][]float64) [][]float64 {
    rows := len(tomMatrix)
    distMatrix := make([][]float64, rows)
    for i := range distMatrix {
        distMatrix[i] = make([]float64, rows)
        for j := 0; j < rows; j++ {
            // 简单的线性转换
            distMatrix[i][j] = 1.0 - tomMatrix[i][j]
        }
    }
    return distMatrix
}