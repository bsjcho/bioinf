package bioinf

import "math"

func Slice2D(rows, cols, defVal int) [][]int {
	m := make([][]int, rows)
	x := make([]int, rows*cols)
	for j := range x {
		x[j] = defVal
	}
	for i := range m {
		m[i], x = x[:cols], x[cols:]
	}
	return m
}

func Max(ints ...int) int {
	max := math.MinInt64
	for _, x := range ints {
		if x > max {
			max = x
		}
	}
	return max
}
