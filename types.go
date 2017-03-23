package bioinf

// Sequence represents a nucleotide base sequence
type Sequence struct {
	Bases []Base
}

// NewSequence is a Sequence constructor
func NewSequence() *Sequence {
	return &Sequence{Bases: []Base{}}
}

// Base represents a nucleotide base
type Base int

// A ... enum represents a nucleotide
const (
	A Base = iota
	C
	G
	T
	X // represents a gap "-"
)

// AsToSeqs converts strings to Sequences
func AsToSeqs(seqStrs []string) (seqs []*Sequence) {
	for _, seqStr := range seqStrs {
		seqs = append(seqs, AToSeq(seqStr))
	}
	return
}

// AToSeq converts string to Sequence
func AToSeq(seq string) *Sequence {
	s := NewSequence()
	for _, b := range seq {
		s.Bases = append(s.Bases, AToBase(string(b)))
	}
	return s
}

// AToBase converts string to Base
func AToBase(b string) Base {
	switch b {
	case "A":
		return A
	case "C":
		return C
	case "G":
		return G
	case "T":
		return T
	default:
		return X
	}
}
