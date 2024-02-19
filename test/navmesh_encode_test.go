package test

import (
	"github.com/gorustyt/gonavmesh/detour"
	"os"
	"testing"
)

func TestNavMeshEncode(t *testing.T) {
	data, err := os.ReadFile("test/test_data/solo_navmesh.bin")
	if err != nil {
		panic(err)
	}
	mesh := &detour.NavMeshData{}
	err = mesh.FromBin(data)
	if err != nil {
		panic(err)
	}
}
