package main

import (
	"bufio"
	"fmt"
	"math"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

type ChrRec struct {
	ChrName         string
	RecName         string
	Seq             string
	SenseEnergy     []float32
	SenseProb       []float32
	AntiSenseEnergy []float32
	AntiSenseProb   []float32
}

type Rec struct {
	Name               string
	Length             int
	EnergyTable        map[string]float32
	RevCompEnergyTable map[string]float32
	ProbTable          map[string]float32
	RevCompProbTable   map[string]float32
}

func (rec *Rec) Scan(chrname, chr string) ChrRec {
	var chrRec ChrRec
	chrRec.ChrName = chrname
	chrRec.RecName = rec.Name
	chrRec.Seq = chr
	chrRec.SenseEnergy = make([]float32, len(chr))
	chrRec.SenseProb = make([]float32, len(chr))
	chrRec.AntiSenseEnergy = make([]float32, len(chr))
	chrRec.AntiSenseProb = make([]float32, len(chr))

	score := make([]float32, len(chr))
	length := 8
	for i := 0; i < (len(chr) - (length - 1)); i++ {
		score[i] = table[chr[i:i+length]]
	}
	// fmt.Println(chr)
	// fmt.Println(score)
	return score
}

func revcomp(s string) string {
	rc := make([]byte, len(s))
	for i, j := 0, len(s)-1; i < len(s); i, j = i+1, j-1 {
		switch s[j] {
		case 'A':
			rc[i] = 'T'
		case 'T':
			rc[i] = 'A'
		case 'C':
			rc[i] = 'G'
		case 'G':
			rc[i] = 'C'
		default:
			rc[i] = s[j]
		}
	}
	return string(rc)
}

func readChr(fn string) string {
	fastaFile, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	var s []alphabet.Letter
	t := linear.NewSeq("", s, alphabet.DNA)
	reader := fasta.NewReader(fastaFile, t)
	seq, _ := reader.Read()
	//fmt.Println("Read -> ", seq.Alphabet())
	return seq.(*linear.Seq).String()
}

func readTable(recname, fn string) *Rec {
	energy := make(map[string]float32)
	prob := make(map[string]float32)
	energyRevComp := make(map[string]float32)
	probRevComp := make(map[string]float32)

	f, err := os.Open(fn)
	if err != nil {
		fmt.Println("ERROR: reading rule", err)
		panic(err)
	}

	scanner := bufio.NewScanner(f)

	var (
		motif        string
		motifRevComp string
		score        float32
		pos          float32
	)
	pos = 1.0
	for scanner.Scan() {
		fmt.Sscanf(scanner.Text(), "%s %f", &motif, &score)
		energy[motif] = score
		prob[motif] = 1.0 - (pos / float32(math.Pow(4.0, float64(len(motif)))))
		motifRevComp = revcomp(motif)
		energyRevComp[motifRevComp] = score
		probRevComp[motifRevComp] = 1.0 - (pos / float32(math.Pow(4.0, float64(len(motif)))))
		pos += 1.0
		// fmt.Println(motif, motifRevComp)
	}

	return &Rec{Name: recname, Length: len(motif), EnergyTable: energy, ProbTable: prob, RevCompEnergyTable: energyRevComp, RevCompProbTable: probRevComp}
}

func main() {
	ch := readChr("/home/jgcarvalho/gocode/src/github.com/jgcarvalho/TFBS-scan/Hs.chrY.fasta")
	receptor := readTable("P10589", "/home/jgcarvalho/gocode/src/github.com/jgcarvalho/TFBS-scan/P10589.csv")
	// fmt.Println(table)
	// scoreEnergy := scan(ch, table)
	chrData := receptor.Scan("chrY", ch)
	// scorePercentil := scan(ch, receptor.ProbTable)
	fmt.Println("Done")
	// fmt.Println(scoreEnergy[100000:100100])
	fmt.Println(scorePercentil[100000:100100])

}
