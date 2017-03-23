package bioinf

var (
	// values doubled to be able to use integers during calculations
	// final result is converted to float then divided by two
	match    = 6
	mismatch = -4
	gap      = -3
)

// SPScore returns the score for sequences
// TODO - handle cases where sequences are of differing lengths
func SPScore(seqs []*Sequence) (score int) {
	for i := range seqs[0].Bases {
		colBases := []Base{}
		for j := range seqs {
			colBases = append(colBases, seqs[j].Bases[i])
		}
		score += ColumnSPScore(colBases)
	}
	return
}

// ColumnSPScore returns the sum-of-pairs score for a column of bases
func ColumnSPScore(bases []Base) (sum int) {
	for i, bi := range bases[:len(bases)-1] {
		for _, bj := range bases[i+1:] {
			sum += PairScore(bi, bj)
		}
	}
	return
}

// PairScore returns the score of a pair of bases (or gap)
func PairScore(b1, b2 Base) int {
	if b1 == X && b2 == X {
		return 0
	}
	if (b1 == X && b2 != X) ||
		(b2 == X && b1 != X) {
		return gap
	}
	if b1 != b2 {
		return mismatch
	}
	return match
}
