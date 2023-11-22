package gui

import (
	"log"
	"os"
	"path/filepath"
)

func GetFileName(dir string) map[string]string {
	ens, err := os.ReadDir(dir)
	if err != nil {
		log.Println(err)
	}
	res := map[string]string{}
	for _, v := range ens {
		if v.IsDir() {
			continue
		}
		res[v.Name()] = filepath.Join(dir, v.Name())
	}
	return res
}
