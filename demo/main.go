package main

import "fmt"

func main() {
	var t [3]float64
	var c = []float64{1, 2, 3}
	copy(t[:], c)
	fmt.Println(t)
}
