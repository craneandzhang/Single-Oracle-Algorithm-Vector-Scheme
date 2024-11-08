package main

import (
	"fmt"
	"math/rand"
	"time"

	"github.com/cloudflare/circl/ecc/bls12381" // Using Cloudflare Circl for elliptic curve operations
)

const N = 1024

var t1, t2 int64

// Initialize pairing parameters
var pairing = bls12381.NewG1()

func main() {
	// Generate random elements
	W := newRandomElement()
	G_vector := generateRandomVector(N)
	a_vector := generateRandomVector(N)
	b_vector := mapToBinaryVector()

	// Compute C
	C := innerProduct(a_vector, b_vector).Mul(W)
	for _, g := range G_vector {
		C = C.Add(g)
	}

	// Measure time for algorithm1
	start := time.Now().UnixMicro()
	result := !algorithm1(G_vector, a_vector, b_vector, C, W, N)
	end := time.Now().UnixMicro()

	fmt.Printf("Algorithm1 total time: %d microseconds\n", end-start)
	t1 = end - start - t2
	fmt.Printf("Algorithm1 execution time: %d microseconds\n", t1)
	fmt.Println("Algorithm1 result:", result)
}

// Main recursive algorithm
func algorithm1(G_vector, a_vector, b_vector []*bls12381.G1, C, W *bls12381.G1, n int) bool {
	if n == 1 {
		start := time.Now().UnixMicro()
		t := innerProduct(a_vector, b_vector)
		tmp := t.Mul(G_vector[0]).Add(a_vector[0].Mul(G_vector[0]))
		result := C.IsEqual(tmp)
		end := time.Now().UnixMicro()
		t2 = end - start
		fmt.Printf("Algorithm1 verification time: %d microseconds\n", t2)
		return result
	}

	mid := n / 2

	// Split arrays
	GL, GR := G_vector[:mid], G_vector[mid:]
	aL, aR := a_vector[:mid], a_vector[mid:]
	bL, bR := b_vector[:mid], b_vector[mid:]

	// Compute L and R
	L := innerProduct(aL, bR).Mul(W).Add(innerProduct(aL, GR))
	R := innerProduct(aR, bL).Mul(W).Add(innerProduct(aR, GL))

	// Generate random t and its inverse
	t := newRandomElement()
	tInverse := newRandomElement().Inverse(t)

	// Update vectors
	a_vector_new := make([]*bls12381.G1, mid)
	b_vector_new := make([]*bls12381.G1, mid)
	G_vector_new := make([]*bls12381.G1, mid)
	for i := 0; i < mid; i++ {
		a_vector_new[i] = aL[i].Add(tInverse.Mul(aR[i]))
		b_vector_new[i] = bL[i].Add(tInverse.Mul(bR[i]))
		G_vector_new[i] = GL[i].Mul(GR[i])
	}

	// Update C for recursion
	C_mew := t.Mul(L).Add(C).Add(tInverse.Mul(R))

	return algorithm1(G_vector_new, a_vector_new, b_vector_new, C_mew, W, mid)
}

// Computes the inner product of two vectors
func innerProduct(a, b []*bls12381.G1) *bls12381.G1 {
	result := pairing.Identity()
	for i := 0; i < len(a); i++ {
		result = result.Add(a[i].Mul(b[i]))
	}
	return result
}

// Generate a random element from Z_p* (non-zero)
func newRandomElement() *bls12381.G1 {
	elem := pairing.Random(rand.Reader)
	for elem.IsIdentity() {
		elem = pairing.Random(rand.Reader)
	}
	return elem
}

// Generates a random vector of elements in Z_p*
func generateRandomVector(size int) []*bls12381.G1 {
	vector := make([]*bls12381.G1, size)
	for i := 0; i < size; i++ {
		vector[i] = newRandomElement()
	}
	return vector
}

// Maps an integer to a binary vector
func mapToBinaryVector() []*bls12381.G1 {
	binaryVector := make([]*bls12381.G1, N)
	for i := 0; i < N; i++ {
		bit := (N >> i) & 1
		if bit == 1 {
			binaryVector[i] = pairing.Identity().Neg()
		} else {
			binaryVector[i] = pairing.Identity()
		}
	}
	return binaryVector
}
