/*
Global pairwise alignment with pair hidden markov models using the viterbi algorithm
Notes:
- matrices in natural log space
- δ=0.08 ε=0.35 τ =0.002
- emission probabilities p and q in their respective .txt files
- proteins being compared are in 2017-01-16uniprot.fasta

To run:
go run viterbi.go

Top 3 alignments comparing first 1000 proteins from 2017-01-16uniprot.fasta
to Z286B_HUMAN:
Index=324 Name=Z286A_HUMAN ln Pr=-2406.2183096044073
METDLAEMPEKGALSSQDSPHFQEKSTEEGEVAALRLTARSQETVTF-----KDVAMDFT
METDLAEMPEKGVLSSQDSPHFQEKSTEEGEVAALRLTARSQAAAAAAAPGSRSLRGVHV
Index=351 Name=ZN419_HUMAN ln Pr=-2875.448973333825
MAAAALRDPAQVPVAADLLTDHEEGYVTFE-DVAVY---F-SQEEWRLLDDAQRLLYR-N
METDLAEMPEK-GVLSSQDSPHFQEKSTEEGEVAALRLTARSQAAAAAAAPGSRSLRGVH
Index=311 Name=ZN157_HUMAN ln Pr=-2880.4033644615038
MPAN-GTSPQRFPALIPGEPGRSF-EGSVSFEDVA-VDFT-R-QEWHRLD-PAQRTM---
METDLAEMPEK--GVLSSQDSPHFQEKSTEEGEVAALRLTARSQAAAAAAAPGSRSLRGV
*/

package main

import (
	"bufio"
	"container/heap"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"time"
)

// Result holds comparison results between a protein and the base
type Result struct {
	index       int
	proteinName string
	lnPrViterbi float64 // natural log of prob. of optimal seq. alignment
	proteinSeq  string
	baseSeq     string
}

// ResultHeap implements heap interface
// ResultHeap is used to maintain top protein matches.
type ResultHeap []*Result

// Protein holds information about a protein to be compared to the base
type Protein struct {
	index int
	name  string
	seq   string
}

// seqState is an enum representing sequence states
type seqState int

const (
	match seqState = iota
	insertion
	deletion
)

const (
	pEmissionsFilename = "p.txt"
	qEmissionsFilename = "q.txt"
	proteinsFilename   = "2017-01-16uniprot.fasta"

	// 1000 Z286B_HUMAN
	base     = `METDLAEMPEKGVLSSQDSPHFQEKSTEEGEVAALRLTARSQAAAAAAAPGSRSLRGVHVPPPLHPAPAREEIKSTCSLKACFSLSLTLTYYRTAFLLSTENEGNLHFQCPSDVETRPQSKDSTSVQDFSKAESCKVAIIDRLTRNSVYDSNLEAALECENWLEKQQGNQERHLREMFTHMNSLSEETDHEHDVYWKSFNQKSVLITEDRVPKGSYAFHTLEKSLKQKSNLMKKQRTYKEKKPHKCNDCGELFTCHSVHIQHQRVHTGEKPYTCNECGKSFSHRANLTKHQRTHTRILFECRECKKTFTESSSLATHQRIHVGERPYECNECGKGFNRSTHLVQHQLIHTGVRPYECNECDKAFIHSSALIKHQRTHTGEKPYKCQECGKAFSHCSSLTKHQRVHTGEKPYECSECGKTFSQSTHLVQHQRIHTGEKPYECSECGKTFSQSSNFAKHQRIHIGKKPYKCSECGKAFIHSSALIQHQRTHTGEKPFRCNECGKSFKCSSSLIRHQRVHTEEQP`
	maxIndex = 999

	// // 59 OXSR1_HUMAN
	// base = `MSEDSSALPWSINRDDYELQEVIGSGATAVVQAAYCAPKKEKVAIKRINLEKCQTSMDELLKEIQAMSQCHHPNIVSYYTSFVVKDELWLVMKLLSGGSVLDIIKHIVAKGEHKSGVLDESTIATILREVLEGLEYLHKNGQIHRDVKAGNILLGEDGSVQIADFGVSAFLATGGDITRNKVRKTFVGTPCWMAPEVMEQVRGYDFKADIWSFGITAIELATGAAPYHKYPPMKVLMLTLQNDPPSLETGVQDKEMLKKYGKSFRKMISLCLQKDPEKRPTAAELLRHKFFQKAKNKEFLQEKTLQRAPTISERAKKVRRVPGSSGRLHKTEDGGWEWSDDEFDEESEEGKAAISQLRSPRVKESISNSELFPTTDPVGTLLQVPEQISAHLPQPAGQIATQPTQVSLPPTAEPAKTAQALSSGSGSQETKIPISLVLRLRNSKKELNDIRFEFTPGRDTAEGVSQELISAGLVDGRDLVIVAANLQKIVEEPQSNRSVTFKLASGVEGSDIPDDGKLIGFAQLSIS`
	// maxIndex = 50

	delta   = 0.08
	epsilon = 0.35
	tau     = 0.002
	// max array sizes for matrices
	maxSeq = 1000
)

var (
	q          map[string]float64
	p          map[string]map[string]float64
	vM, vX, vY [maxSeq][maxSeq]float64 // viterbi matrices
	pro        *Protein                // protein to be compared to base
	results    *ResultHeap
	nINF       = math.Inf(-1)
)

func init() {
	results = &ResultHeap{}
	heap.Init(results)
}

// Program entry point
func main() {
	q = parseQ(qEmissionsFilename)
	p = parseP(pEmissionsFilename)
	initializeTables() // viterbi initialization
	proteins := parseProteins(proteinsFilename)
	start := time.Now()
	for i, protein := range proteins {
		if i > maxIndex {
			break
		}
		pro = protein
		compareProtein()
	}
	duration := time.Since(start)
	fmt.Printf("Execution time: %vs\n", duration.Seconds())
	fmt.Println("Top 3 Results:")
	for i := 0; i < 3; i++ {
		r := heap.Pop(results)
		result, _ := r.(*Result)
		result.print()
	}
}

// viterbi initialization
// only runs once
func initializeTables() {
	vM[1][1] = math.Log(1)
	vX[1][1] = nINF
	vY[1][1] = nINF
	for i := 0; i < maxSeq; i++ {
		updateTables(i, 0, nINF)
		updateTables(0, i, nINF)
	}
}

// implicit comparison to global var pro :\
func compareProtein() {
	// fill viterbi DP tables (vm, vx, vy)
	fillTables()
	result := traceback()
	result.index = pro.index
	result.proteinName = pro.name
	heap.Push(results, result) // using heap to maintain order of results
	result.print()
}

// reference materials use negative indices (specifically, -1) in matrices.
// arrays with negative indices are not permissible in golang.
// so i've incremented all indices by 1.
// -1 is now represented as 0
// adjustments have been made throughout the code to account for this :(
// TODO: produce an accessor for the matrices to handle negative index access
func fillTables() {
	for i := 1; i <= lenProSeq(); i++ {
		for j := 1; j <= lenBaseSeq(); j++ {
			if i == 1 && j == 1 {
				continue
			}
			// viterbi recurrence
			vM[i][j] = getP(i, j) + math.Max(
				vM[i-1][j-1]+math.Log(1-(2*delta)-tau),
				// if only golang had a built-in variadic Max function :\
				math.Max(
					vX[i-1][j-1]+math.Log(1-epsilon-tau),
					vY[i-1][j-1]+math.Log(1-epsilon-tau),
				),
			)
			vX[i][j] = getQ(i, pro.seq) + math.Max(
				vM[i-1][j]+math.Log(delta),
				vX[i-1][j]+math.Log(epsilon),
			)
			vY[i][j] = getQ(j, base) + math.Max(
				vM[i][j-1]+math.Log(delta),
				vY[i][j-1]+math.Log(epsilon),
			)
		}
	}
}

func traceback() *Result {
	result := &Result{}
	max, state := termination()
	result.lnPrViterbi = math.Log(tau) + max
	stateSequence := tracebackHelper(state, len(pro.seq)+1, len(base)+1)
	generateSequences(result, stateSequence)
	return result
}

func termination() (max float64, state seqState) {
	n, m := lenProSeq(), lenBaseSeq()
	vm := vM[n][m]
	vx := vX[n][m]
	vy := vY[n][m]
	max, state = getMaxAndState(vm, vx, vy)
	return
}

// Recursive helper function for traceback
func tracebackHelper(state seqState, i, j int) []seqState {
	states := []seqState{}
	// base case
	if i == 1 && j == 1 {
		return states
	}
	states = append(states, state)
	priorState, pi, pj := determinePriorState(state, i, j)
	return merge(tracebackHelper(priorState, pi, pj), states)
}

// helper function used by recursive traceback function to determine
// the next state and indices to go to
func determinePriorState(state seqState, i, j int) (priorState seqState,
	pi int, pj int) {
	vmVal, vxVal, vyVal := nINF, nINF, nINF
	switch state {
	case match:
		vmVal = vM[i-1][j-1]
		vxVal = vX[i-1][j-1]
		vyVal = vY[i-1][j-1]
		pi = i - 1
		pj = j - 1
	case insertion:
		vmVal = vM[i-1][j]
		vxVal = vX[i-1][j]
		pi = i - 1
		pj = j
	case deletion:
		vmVal = vM[i][j-1]
		vyVal = vY[i][j-1]
		pi = i
		pj = j - 1
	}
	_, priorState = getMaxAndState(vmVal, vxVal, vyVal)
	return
}

// helper function to extract max value and corresponding sequence state
func getMaxAndState(vm, vx, vy float64) (max float64, state seqState) {
	max = nINF
	if vm > max {
		max = vm
		state = match
	}
	if vx > max {
		max = vx
		state = insertion
	}
	if vy > max {
		max = vy
		state = deletion
	}
	return
}

// Generates string representations of optimal alignment
func generateSequences(result *Result, stateSequence []seqState) {
	var proSeq, baseSeq string
	var i, j int
	for _, state := range stateSequence {
		switch state {
		case match:
			proSeq += string(pro.seq[i])
			baseSeq += string(base[j])
			i++
			j++
		case insertion:
			proSeq += string(pro.seq[i])
			baseSeq += "-"
			i++
		case deletion:
			proSeq += "-"
			baseSeq += string(base[j])
			j++
		}
	}
	result.proteinSeq = proSeq
	result.baseSeq = baseSeq
}

///////////////// Helper Functions

func lenProSeq() int {
	return len(pro.seq) + 1
}

func lenBaseSeq() int {
	return len(base) + 1
}

func getP(x, y int) float64 {
	// necessary to decrement input values by 2 due to indexing of viterbi DP tables
	x -= 2
	y -= 2
	if x < 0 || y < 0 {
		return nINF
	}
	return p[string(pro.seq[x])][string(base[y])]
}

func getQ(i int, seq string) float64 {
	// necessary to decrement input values by 2 due to indexing of viterbi DP tables
	i -= 2
	if i < 0 {
		return nINF
	}
	return q[string(seq[i])]
}

func updateTables(i, j int, val float64) {
	vM[i][j] = val
	vX[i][j] = val
	vY[i][j] = val
}

func merge(states, statesToBeMerged []seqState) []seqState {
	for _, state := range statesToBeMerged {
		states = append(states, state)
	}
	return states
}

func convertToLnFloat(floatString string) float64 {
	prob, err := strconv.ParseFloat(floatString, 64)
	checkErr(err)
	return math.Log(prob)
}

///////////////// Parse Functions

func parseProteins(filename string) (proteins []*Protein) {
	scanner, file := fetchScanner(filename)
	defer file.Close()
	var protein *Protein
	index := 0
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		if strings.Contains(fields[0], ">") {
			if index > 0 {
				proteins = append(proteins, protein)
			}
			protein = &Protein{}
			protein.index = index
			protein.name = strings.Split(fields[0], "|")[2]
			index++
		} else {
			protein.seq += strings.Replace(strings.TrimSpace(fields[0]), "U", "C", -1)
		}
	}
	return
}

func parseP(filename string) map[string]map[string]float64 {
	rp := map[string]map[string]float64{}
	scanner, file := fetchScanner(filename)
	defer file.Close()
	var order []string
	firstLine := true
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		if firstLine {
			order = fields
			firstLine = false
		} else {
			var xKey, yKey string
			py := map[string]float64{}
			for i, val := range fields {
				if i == 0 {
					xKey = val
				} else {
					yKey = order[i-1]
					py[yKey] = convertToLnFloat(val)
				}
			}
			rp[xKey] = py
		}
	}
	return rp
}

func parseQ(filename string) map[string]float64 {
	rq := map[string]float64{}
	scanner, file := fetchScanner(filename)
	defer file.Close()
	for scanner.Scan() {
		s := strings.Fields(scanner.Text())
		rq[s[0]] = convertToLnFloat(s[1])
	}
	return rq
}

func fetchScanner(filename string) (*bufio.Scanner, *os.File) {
	file, err := os.Open(filename)
	checkErr(err)
	return bufio.NewScanner(file), file
}

////////////////// ResultsHeap heap interface implementation

func (h ResultHeap) Len() int {
	return len(h)
}

func (h ResultHeap) Less(i, j int) bool {
	return h[i].lnPrViterbi > h[j].lnPrViterbi
}

func (h ResultHeap) Swap(i, j int) {
	h[i], h[j] = h[j], h[i]
}

// Push - heap interface
func (h *ResultHeap) Push(x interface{}) {
	*h = append(*h, x.(*Result))
}

// Pop - heap interface
func (h *ResultHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

////////////////// Aux Functions

func (r *Result) print() {
	fmt.Printf("Index=%v Name=%v ln Pr=%v\n%v\n%v\n", r.index,
		r.proteinName, r.lnPrViterbi, r.proteinSeq[:60], r.baseSeq[:60])
}

func checkErr(err error) {
	if err != nil {
		fmt.Printf("***** ERROR *****\n%v\n", err)
	}
}
