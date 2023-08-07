// Copyright (c) 2023 Colin McRae

package bigmatrix

import (
	"fmt"
	"github.com/stretchr/testify/assert"
	"os"
	"pslq/bignumber"
	"testing"
)

func TestMain(m *testing.M) {
	err := bignumber.Init(500)
	if err != nil {
		fmt.Printf("Invalid input to Init: %q", err.Error())
		return
	}
	code := m.Run()
	os.Exit(code)
}

func TestNewFromDecimalStringArray(t *testing.T) {
	x, err := NewFromDecimalStringArray([]string{"0"}, 1, 2)
	assert.Error(t, err)
	assert.Nil(t, x)

	x, err = NewFromDecimalStringArray([]string{}, 0, 1)
	assert.Error(t, err)
	assert.Nil(t, x)

	x, err = NewFromDecimalStringArray([]string{"0", "0", "a"}, 3, 1)
	assert.Error(t, err)
	assert.Nil(t, x)
}

func TestNewIdentity(t *testing.T) {
	identity, err := NewIdentity(3)
	assert.NoError(t, err)
	assert.NotNil(t, identity)
	assert.Equal(t, 3, identity.numRows)
	assert.Equal(t, 3, identity.numCols)
	zero := bignumber.NewFromInt64(0)
	one := bignumber.NewFromInt64(1)
	for i := 0; i < identity.numRows; i++ {
		for j := 0; j < identity.numCols; j++ {
			if i == j {
				assert.Equal(t, 0, identity.values[i*identity.numCols+j].Cmp(one))
			} else {
				assert.Equal(t, 0, identity.values[i*identity.numCols+j].Cmp(zero))
			}
		}
	}

	// Dimension 0 or less
	_, err = NewIdentity(0)
	assert.Error(t, err)
}

func TestBigMatrix_Add(t *testing.T) {
	testBinaryOp(t, //          1           2          3           4           5           6         7            8         9           10            11       12          13      14          15          16         17         18        19
		[]string{"30.9526", "-23.7367", "-13.4514", "-2.2083", "-14.642", "-20.4062", "-43.2513", "82.8234", "-36.4775", "17.094", "9.9718", "-41.9313", "-22.9266", "74.152", "-23.184", "-22.5942", "107.7407", "-70.3026", "7.3085", "11.62"},
		4, 5, //      1             2           3            4      5         6               7        8          9       10               11        12          13         14          15          16        17          18           19
		[]string{"-1572.98", "-165884.2623", "-1345.14", "-4827.59", "198", "3708.12", "16228.8163", "8282.34", "2282.21", "-8239", "49471.04", "-2497.11", "-229266", "12527.686", "-231.84", "151316.7", "-7487.27", "-703026", "38256.9175", "-116.2"},
		4, 5, //           1           2            3            4            5            6               7           8         9           10              11              12          13           14            15            16              17            18           19
		[]string{"-1542.0274", "-165907.999", "-1358.5914", "-4829.7983", "183.358", "3687.7138", "16185.565", "8365.1634", "2245.7325", "-8221.906", "49481.0118", "-2539.0413", "-229288.9266", "12601.838", "-255.024", "151294.1058", "-7379.5293", "-703096.3026", "38264.226", "-104.58"},
		4, 5,
		"Add",
	)
}

func TestBigMatrix_Sub(t *testing.T) {
	testBinaryOp(t, //          1           2          3           4           5           6         7            8         9           10            11       12          13      14          15          16         17          18       19
		[]string{"30.9526", "-23.7367", "-13.4514", "-2.2083", "-14.642", "-20.4062", "-43.2513", "82.8234", "-36.4775", "17.094", "9.9718", "-41.9313", "-22.9266", "74.152", "-23.184", "-22.5942", "107.7407", "-70.3026", "7.3085", "11.62"},
		4, 5, //      1             2           3            4      5         6               7        8          9       10               11        12          13         14          15          16         17         18             19
		[]string{"-1572.98", "-165884.2623", "-1345.14", "-4827.59", "198", "3708.12", "16228.8163", "8282.34", "2282.21", "-8239", "49471.04", "-2497.11", "-229266", "12527.686", "-231.84", "151316.7", "-7487.27", "-703026", "38256.9175", "-116.2"},
		4, 5, //           1           2            3            4            5            6               7           8         9           10              11              12          13           14            15            16              17             18         19
		[]string{"1603.9326", "165860.5256", "1331.6886", "4825.3817", "-212.642", "-3728.5262", "-16272.0676", "-8199.5166", "-2318.6875", "8256.094", "-49461.0682", "2455.1787", "229243.0734", "-12453.534", "208.656", "-151339.2942", "7595.0107", "702955.6974", "-38249.609", "127.82"},
		4, 5,
		"Sub",
	)
}

func TestBigMatrix_Mul(t *testing.T) {
	testBinaryOp(t,
		[]string{"3.14", "-1.59", "2.65", "-3.5", "9.79", "-3.23", "0", "-2.71", "8.28", "0", "-8.31", "-4.15"},
		4, 3,
		[]string{"9.26", "-5.35", "0", "-9.32", "-2.3", "1.82", "-8.97", "8.46", "-4.6", "0", "1.8", "-8", "0", "7.45", "-2.8"},
		3, 5,
		[]string{"30.9526", "-23.7367", "-13.4514", "-2.2083", "-14.642", "-20.4062", "-43.2513", "82.8234", "-36.4775", "17.094", "9.9718", "-41.9313", "-22.9266", "74.152", "-23.184", "-22.5942", "107.7407", "-70.3026", "7.3085", "11.62"},
		4, 5,
		"Mul",
	)
	testBinaryOp(t,
		[]string{"-314", "1.59", "265", "35", "-9.79", "-323", "0", "271", "-8.28", "0", "831", "-4.15"},
		4, 3,
		[]string{"-9.26", "535", "0", "9.32", "23", "182", "-8.97", "-846", "46", "0", "-18", "8", "0", "-7.45", "28"},
		3, 5,
		[]string{"-1572.98", "-165884.2623", "-1345.14", "-4827.59", "198", "3708.12", "16228.8163", "8282.34", "2282.21", "-8239", "49471.04", "-2497.11", "-229266", "12527.686", "-231.84", "151316.7", "-7487.27", "-703026", "38256.9175", "-116.2"},
		4, 5,
		"Mul",
	)

	// Create matrices needed below
	oneByTwo, err := NewFromDecimalStringArray([]string{"1", "2"}, 1, 2)
	assert.NoError(t, err)
	wrongLen, err := NewFromDecimalStringArray([]string{"1", "2", "3", "4"}, 2, 2)
	assert.NoError(t, err)
	wrongLen.numRows = 1
	negativeDim, err := NewFromDecimalStringArray([]string{"1", "2", "3", "4"}, 2, 2)
	negativeDim.numRows = -2
	negativeDim.numCols = -2
	assert.NoError(t, err)

	// Mismatched dimensions
	out, err := NewEmpty(0, 0).Mul(oneByTwo, oneByTwo)
	if out != nil {
		assert.True(t, false)
	}
	assert.Error(t, err)

	// In x, length of input != product of dimensions
	out, err = NewEmpty(0, 0).Mul(wrongLen, oneByTwo)
	if out != nil {
		assert.True(t, false)
	}
	assert.Error(t, err)

	// In y, negative dimensions
	out, err = NewEmpty(0, 0).Mul(oneByTwo, negativeDim)
	if out != nil {
		assert.True(t, false)
	}
	assert.Error(t, err)

	// In y, length of input != product of dimensions
	out, err = NewEmpty(0, 0).Mul(oneByTwo, wrongLen)
	if out != nil {
		assert.True(t, false)
	}
	assert.Error(t, err)

	// In x, negative dimensions
	out, err = NewEmpty(0, 0).Mul(negativeDim, oneByTwo)
	if out != nil {
		assert.True(t, false)
	}
	assert.Error(t, err)
}

func TestBigMatrix_Int64DotProduct(t *testing.T) {
	const (
		xNumRows = 4
		xNumCols = 3
		yNumCols = 5
	)
	xEntries := []int64{
		1, -3, 5,
		-7, 2, -4,
		8, -16, 13,
		-12, 8, 9,
	}
	yEntries := []string{
		"3.14", "2.71", "-82.8", "-0.159", "2.65",
		"-35", "8.9", "-18.2", "-50.2", "88.9",
		"-27.1", "44.6", "-91", "-.0123", "1.333",
	}
	expectedXYEntries := []string{
		"-27.36", "199.01", "-483.2", "150.3795", "-257.385",
		"16.42", "-179.57", "907.2", "-99.2378", "153.918",
		"232.82", "459.08", "-1554.2", "801.7681", "-1383.871",
		"-561.58", "440.08", "29", "-399.8027", "691.397",
	}
	y, err := NewFromDecimalStringArray(yEntries, xNumCols, yNumCols)
	assert.NoError(t, err)
	actualXY := NewEmpty(xNumRows, yNumCols)
	expectedXY, err := NewFromDecimalStringArray(expectedXYEntries, xNumRows, yNumCols)
	assert.NoError(t, err)
	for i := 0; i < xNumRows; i++ {
		for j := 0; j < yNumCols; j++ {
			entry, err := Int64DotProduct(xEntries, xNumCols, y, i, j, 0, xNumCols, false)
			assert.NoError(t, err)
			err = actualXY.Set(i, j, entry)
			assert.NoError(t, err)
		}
	}
	equals, err := expectedXY.Equals(actualXY, bignumber.NewPowerOfTwo(-50))
	assert.True(t, equals)
}

func TestBigMatrix_Equals(t *testing.T) {
	var errMessage string
	log2tolerance := -50
	tolerance := bignumber.NewPowerOfTwo(int64(log2tolerance))
	nearDelta := bignumber.NewPowerOfTwo(int64(log2tolerance - 1))
	mediumDelta := bignumber.NewPowerOfTwo(int64(log2tolerance + 1))
	s := []string{
		"1.23", "4.56", "0.12", "5.67",
		"12.3", "45.6", "10.12", "51.67",
		"111.235", "45.6", "41.12", "9.67",
	}
	a, err := NewFromDecimalStringArray(s, 3, 4)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "could not create a: %q", errMessage)

	// Test close but within tolerance
	aNear, err := NewFromDecimalStringArray(s, 3, 4)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "could not create aNear: %q", errMessage)
	aNear.values[5].Add(aNear.values[5], nearDelta)
	aNear.values[7].Add(aNear.values[7], nearDelta)
	aNear.values[9].Add(aNear.values[9], nearDelta)
	aNear.values[11].Add(aNear.values[11], nearDelta)
	equals, err := a.Equals(aNear, tolerance)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error in Equals: %q", errMessage)
	assert.True(t, equals)

	// Test just outside of tolerance
	aMedium, err := NewFromDecimalStringArray(s, 3, 4)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "could not create aNear: %q", errMessage)
	aMedium.values[3].Add(aNear.values[3], mediumDelta)
	aMedium.values[6].Add(aNear.values[6], nearDelta)
	aMedium.values[8].Add(aNear.values[8], nearDelta)
	aMedium.values[10].Add(aNear.values[10], nearDelta)
	equals, err = a.Equals(aMedium, tolerance)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error in Equals: %q", errMessage)
	assert.False(t, equals)

	// Mismatched dimensions
	twoByThree, err := NewFromDecimalStringArray([]string{"1", "2", "3", "4", "5", "6"}, 2, 3)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "could not create twoByThree: %q", errMessage)
	twoByTwo, err := NewFromDecimalStringArray([]string{"1", "2", "3", "4"}, 2, 2)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "could not create twoByThree: %q", errMessage)
	equals, err = twoByThree.Equals(twoByTwo, bignumber.NewPowerOfTwo(0))
	assert.False(t, equals)
	assert.Error(t, err)

	// Empty input
	empty := NewEmpty(0, 0)
	equals, err = empty.Equals(NewEmpty(0, 0), bignumber.NewFromInt64(0))
	assert.NoError(t, err)
	assert.True(t, equals)

	// Empty and non-empty input
	equals, err = empty.Equals(twoByThree, bignumber.NewFromInt64(0))
	assert.Error(t, err)
	assert.False(t, equals)
}

func TestBigMatrix_Dimensions(t *testing.T) {
	x := NewEmpty(4, 5)
	numRows, numCols := x.Dimensions()
	assert.Equal(t, 4, numRows)
	assert.Equal(t, 5, numCols)
}

func TestBigMatrix_Set(t *testing.T) {
	actual := NewEmpty(3, 5)
	for i := 0; i < 3; i++ {
		for j := 0; j < 5; j++ {
			// 0, 3, 6, 9, 12, 15, ..., 42
			err := actual.Set(i, j, bignumber.NewFromInt64(int64(15*i+3*j)))
			assert.NoError(t, err)
		}
	}
	expected, err := NewFromInt64Array(
		[]int64{0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42}, 3, 5,
	)
	assert.NoError(t, err)
	equals, err := expected.Equals(actual, bignumber.NewFromInt64(0))
	assert.NoError(t, err)
	assert.True(t, equals)

	// Row index out of range
	err = actual.Set(3, 0, bignumber.NewFromInt64(0))
	assert.Error(t, err)
	err = actual.Set(-1, 0, bignumber.NewFromInt64(0))
	assert.Error(t, err)

	// Column index out of range
	err = actual.Set(0, 5, bignumber.NewFromInt64(0))
	assert.Error(t, err)
	err = actual.Set(0, -1, bignumber.NewFromInt64(0))
	assert.Error(t, err)
}

func TestBigMatrix_Get(t *testing.T) {
	x := NewEmpty(3, 5)
	for i := 0; i < 3; i++ {
		for j := 0; j < 5; j++ {
			// 0, 3, 6, 9, 12, 15, ..., 42
			err := x.Set(i, j, bignumber.NewFromInt64(int64(15*i+3*j)))
			assert.NoError(t, err)
		}
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 5; j++ {
			actual, err := x.Get(i, j)
			assert.NoError(t, err)
			assert.Equal(t, 0, actual.Cmp(bignumber.NewFromInt64(int64(15*i+3*j))))
		}
	}

	// Row index out of range
	y1, err := x.Get(3, 0)
	assert.Nil(t, y1)
	assert.Error(t, err)
	y2, err := x.Get(-1, 0)
	assert.Nil(t, y2)
	assert.Error(t, err)

	// Column index out of range
	y3, err := x.Get(0, 5)
	assert.Nil(t, y3)
	assert.Error(t, err)
	y4, err := x.Get(0, -1)
	assert.Nil(t, y4)
	assert.Error(t, err)
}

func TestBigMatrix_Transpose(t *testing.T) {
	// 1 x 4 and 4 x 1
	testTranspose(
		t,
		[]int64{1, 3, 5, 7}, 1, 4,
		[]int64{1, 3, 5, 7}, 4, 1,
	)
	testTranspose(
		t,
		[]int64{1, 3, 5, 7}, 4, 1,
		[]int64{1, 3, 5, 7}, 1, 4,
	)

	// 5x3 and 3x5
	//  1 3 5
	//  7 9 2       1 7 4 10 7
	//  4 6 8   ->  3 9 6  8 5
	// 10 8 7       5 2 8  7 3
	//  7 5 3
	testTranspose(
		t,
		[]int64{1, 3, 5, 7, 9, 2, 4, 6, 8, 10, 8, 7, 7, 5, 3}, 5, 3,
		[]int64{1, 7, 4, 10, 7, 3, 9, 6, 8, 5, 5, 2, 8, 7, 3}, 3, 5,
	)
	testTranspose(
		t,
		[]int64{1, 7, 4, 10, 7, 3, 9, 6, 8, 5, 5, 2, 8, 7, 3}, 3, 5,
		[]int64{1, 3, 5, 7, 9, 2, 4, 6, 8, 10, 8, 7, 7, 5, 3}, 5, 3,
	)

	// 4 x 4
	//  1  3  5  7
	//  2  4  6  8
	//  4  9 16 25
	// -1 -2 -3 -4
	testTranspose(
		t,
		[]int64{1, 3, 5, 7, 2, 4, 6, 8, 4, 9, 16, 25, -1, -2, -3, -4}, 4, 4,
		[]int64{1, 2, 4, -1, 3, 4, 9, -2, 5, 6, 16, -3, 7, 8, 25, -4}, 4, 4,
	)
}

func testBinaryOp(
	t *testing.T,
	xStrs []string,
	xRows int,
	xCols int,
	yStrs []string,
	yRows int,
	yCols int,
	expectedStrs []string,
	expectedRows int,
	expectedCols int,
	operation string,
) {
	var errMessage string
	x, err := NewFromDecimalStringArray(xStrs, xRows, xCols)
	if err != nil {
		errMessage = err.Error()
	}
	xAsReceiver, err := NewFromDecimalStringArray(xStrs, xRows, xCols)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error when creating x: %q", errMessage)
	y, err := NewFromDecimalStringArray(yStrs, yRows, yCols)
	if err != nil {
		errMessage = err.Error()
	}
	yAsReceiver, err := NewFromDecimalStringArray(yStrs, yRows, yCols)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error when creating y: %q", errMessage)
	expected, err := NewFromDecimalStringArray(expectedStrs, expectedRows, expectedCols)
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error when creating expected: %q", errMessage)
	var actual0, actual1, actual2 *BigMatrix

	if operation == "Mul" {
		actual0, err = NewEmpty(0, 0).Mul(x, y)
		actual1, err = xAsReceiver.Mul(xAsReceiver, y)
		actual2, err = yAsReceiver.Mul(x, yAsReceiver)
	}
	if operation == "Add" {
		actual0, err = NewEmpty(0, 0).Add(x, y)
		actual1, err = xAsReceiver.Add(xAsReceiver, y)
		actual2, err = yAsReceiver.Add(x, yAsReceiver)
	}
	if operation == "Sub" {
		actual0, err = NewEmpty(0, 0).Sub(x, y)
		actual1, err = xAsReceiver.Sub(xAsReceiver, y)
		actual2, err = yAsReceiver.Sub(x, yAsReceiver)
	}
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from %s: %q", operation, errMessage)

	// Test actual0
	checkBinaryOutput(t, expected, actual0)
	checkBinaryOutput(t, expected, actual1)
	checkBinaryOutput(t, expected, actual2)
	checkBinaryOutput(t, expected, xAsReceiver)
	checkBinaryOutput(t, expected, yAsReceiver)
}

func checkBinaryOutput(t *testing.T, expected, actual *BigMatrix) {
	var errMessage string

	equals, err := expected.Equals(actual, bignumber.NewPowerOfTwo(-50))
	if err != nil {
		errMessage = err.Error()
	}
	assert.NoErrorf(t, err, "unexpected error from Equals: %q", errMessage)
	assert.True(t, equals)

}

func testTranspose(
	t *testing.T,
	inputAsInt64 []int64, inputRows int, inputCols int,
	expectedAsInt64 []int64, expectedRows int, expectedCols int,
) {
	inputAsMatrix, err := NewFromInt64Array(inputAsInt64, inputRows, inputCols)
	assert.NoError(t, err)
	expectedAsMatrix, err := NewFromInt64Array(expectedAsInt64, expectedRows, expectedCols)
	assert.NoError(t, err)

	// Using any matrix other than the input as the receiver should work
	actualA, err := NewEmpty(5, 3).Transpose(inputAsMatrix) // initial dimensions don't matter
	assert.NoError(t, err)
	equals, err := expectedAsMatrix.Equals(actualA, bignumber.NewFromInt64(0))
	assert.NoError(t, err)

	// Using input as the receiver should convert the input to its transpose
	actualB, err := inputAsMatrix.Transpose(inputAsMatrix)
	assert.NoError(t, err)
	equals, err = expectedAsMatrix.Equals(actualB, bignumber.NewFromInt64(0))
	assert.NoError(t, err)
	equals, err = expectedAsMatrix.Equals(inputAsMatrix, bignumber.NewFromInt64(0))
	assert.NoError(t, err)
	assert.True(t, equals)
}
