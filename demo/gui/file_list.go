package gui

import (
	"log"
	"os"
	"path/filepath"
)

func GetFileName(dir string, ext string, out map[string]string) {
	ens, err := os.ReadDir(dir)
	if err != nil {
		log.Println(err)
	}
	for _, v := range ens {
		if v.IsDir() || filepath.Ext(v.Name()) != ext {
			continue
		}
		out[v.Name()] = filepath.Join(dir, v.Name())
	}
}
