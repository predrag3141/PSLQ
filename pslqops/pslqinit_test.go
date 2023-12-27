// Copyright (c) 2023 Colin McRae

package pslqops

import (
	"testing"

	"github.com/predrag3141/PSLQ/bigmatrix"
	"github.com/predrag3141/PSLQ/bignumber"
	"github.com/stretchr/testify/assert"
)

func TestGetNormalizedX(t *testing.T) {
	// Known input and output from pslq.py
	input := []string{"1", "-3", "-5", "70"}
	expected, err := bigmatrix.NewFromDecimalStringArray(
		[]string{
			"0.014234965584343246", "-0.04270489675302974",
			"-0.07117482792171623", "0.9964475909040272",
		}, 1, 4,
	)
	assert.NoError(t, err)
	rawX, err := GetRawX(input)
	assert.NoError(t, err)
	actual, err := GetNormalizedX(rawX)
	assert.NoError(t, err)
	equals, err := expected.Equals(actual, bignumber.NewPowerOfTwo(-50))
	assert.True(t, equals)

	// Input cannot be parsed
	input = []string{"0.23402694307324226", "a", "0.7352167040894672"}
	shouldBeNil, err := GetRawX(input)
	assert.Error(t, err)
	assert.Nil(t, shouldBeNil)

	// Input contains a zero
	input = []string{"0.23402694307324226", "0", "0.7352167040894672"}
	shouldBeNil, err = GetRawX(input)
	assert.Error(t, err)
	assert.Nil(t, shouldBeNil)
}

func TestGetS(t *testing.T) {
	// Known input and output from pslq.py
	input := []string{"0.23402694307324226", "0.6361507588171331", "0.7352167040894672"}
	expected, err := bigmatrix.NewFromDecimalStringArray(
		[]string{"1", "0.972230111607223", "0.7352167040894672"}, 1, 3,
	)
	assert.NoError(t, err)
	actual, err := GetSFromStrs(input)
	assert.NoError(t, err)
	equals, err := expected.Equals(actual, bignumber.NewPowerOfTwo(-50))
	assert.True(t, equals)

	// Input cannot be parsed
	input = []string{"0.23402694307324226", "a", "0.7352167040894672"}
	shouldBeNil, err := GetSFromStrs(input)
	assert.Error(t, err)
	assert.Nil(t, shouldBeNil)

	// Input contains a zero
	input = []string{"0.23402694307324226", "0", "0.7352167040894672"}
	shouldBeNil, err = GetSFromStrs(input)
	assert.Error(t, err)
	assert.Nil(t, shouldBeNil)

	// Input is not a row vector
	shouldBeNil, err = GetS(bigmatrix.NewEmpty(3, 1))
	assert.Error(t, err)
	assert.Nil(t, shouldBeNil)
}

func TestGetH(t *testing.T) {
	// Known input and output from pslq.py
	normalizedX0, err := bigmatrix.NewFromDecimalStringArray(
		[]string{"0.23402694307324226", "0.6361507588171331", "0.7352167040894672"},
		1, 3,
	)
	assert.NoError(t, err)
	s0, err := bigmatrix.NewFromDecimalStringArray(
		[]string{"1", "0.972230111607223", "0.7352167040894672"}, 1, 3,
	)
	assert.NoError(t, err)
	expectedH0, err := bigmatrix.NewFromDecimalStringArray(
		[]string{
			"0.9722301116072231", "0",
			"-0.15312878673710792", "0.7562167590901482",
			"-0.17697509643062187", "-0.6543211851004007",
		},
		3, 2,
	)
	assert.NoError(t, err)
	actualH0, err := GetH(normalizedX0, s0)
	assert.NoError(t, err)
	equals, err := expectedH0.Equals(actualH0, bignumber.NewPowerOfTwo(-50))
	assert.True(t, equals)
	testPropertiesOfH(t, normalizedX0, actualH0)

	// Generate a new H and test its properties
	rawX1, err := GetRawX([]string{"1", "2435", "3524", "23452", "674"})
	assert.NoError(t, err)
	normalizedX1, err := GetNormalizedX(rawX1)
	assert.NoError(t, err)
	s1, err := GetS(normalizedX1)
	assert.NoError(t, err)
	actualH1, err := GetH(normalizedX1, s1)
	testPropertiesOfH(t, rawX1, actualH1)
}

func testPropertiesOfH(t *testing.T, x *bigmatrix.BigMatrix, h *bigmatrix.BigMatrix) {
	// H-transpose H = I(n-1)
	ht, err := bigmatrix.NewEmpty(0, 0).Transpose(h)
	assert.NoError(t, err)
	shouldBeIdentity, err := bigmatrix.NewEmpty(0, 0).Mul(ht, h)
	assert.NoError(t, err)
	isTheIdentity, err := bigmatrix.NewIdentity(h.NumCols())
	assert.NoError(t, err)
	equals, err := isTheIdentity.Equals(shouldBeIdentity, bignumber.NewPowerOfTwo(-50))
	assert.True(t, equals)

	// x H = 0
	zero := bigmatrix.NewEmpty(1, h.NumCols())
	shouldBeZero, err := bigmatrix.NewEmpty(0, 0).Mul(x, h)
	assert.NoError(t, err)
	equals, err = zero.Equals(shouldBeZero, bignumber.NewPowerOfTwo(-50))
	assert.NoError(t, err)
	assert.True(t, equals)
}
