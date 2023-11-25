package common

func AssertTrue(value bool) {
	if !value {
		panic("Assert not true")
	}
}
