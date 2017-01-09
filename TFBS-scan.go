package main

import (
	"bufio"
	"flag"
	"fmt"
	"math"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

// Estrutura que armazena os dados do cromosso e os valores de energia da ligação com o receptor para cada base. O valor da ligação é anotado na primeira base do motivo
type ChrRec struct {
	ChrName         string
	RecName         string
	RecLength       int
	Seq             string
	SenseEnergy     []float32
	SenseProb       []float32
	AntiSenseEnergy []float32
	AntiSenseProb   []float32
}

// Estrutura do receptor. Nome, comprimento e tabelas de energia e probabilidade para a fita senso e antisenso
type Rec struct {
	Name               string
	Length             int
	EnergyTable        map[string]float32
	RevCompEnergyTable map[string]float32
	ProbTable          map[string]float32
	RevCompProbTable   map[string]float32
}

// type Peaks struct {
// }

// Função de scan que atribui os valores de energia e probabilidade para cada base do cromossomo nas fitas senso e antisenso
func (rec *Rec) Scan(chrName, chr string) ChrRec {
	var chrRec ChrRec
	chrRec.ChrName = chrName
	chrRec.RecName = rec.Name
	chrRec.RecLength = rec.Length
	chrRec.Seq = chr
	chrRec.SenseEnergy = make([]float32, len(chr))
	chrRec.SenseProb = make([]float32, len(chr))
	chrRec.AntiSenseEnergy = make([]float32, len(chr))
	chrRec.AntiSenseProb = make([]float32, len(chr))

	length := rec.Length
	for i := 0; i < (len(chr) - (length - 1)); i++ {
		chrRec.SenseEnergy[i] = rec.EnergyTable[chr[i:i+length]]
		chrRec.SenseProb[i] = rec.ProbTable[chr[i:i+length]]
		chrRec.AntiSenseEnergy[i] = rec.RevCompEnergyTable[chr[i:i+length]]
		chrRec.AntiSenseProb[i] = rec.RevCompProbTable[chr[i:i+length]]
	}

	return chrRec
}

// complemento reverso
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

// Le o fasta do cromossomo
func readChr(fn string) (string, string) {
	fastaFile, err := os.Open(fn)
	if err != nil {
		fmt.Println("Erro ao ler o arquivo", err)
	}
	var s []alphabet.Letter
	t := linear.NewSeq("", s, alphabet.DNA)
	reader := fasta.NewReader(fastaFile, t)
	seq, _ := reader.Read()
	return seq.Name(), seq.(*linear.Seq).String()
}

// Le a tabela de energias ordenada em ordem decrescente
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
	}

	return &Rec{Name: recname, Length: len(motif), EnergyTable: energy, ProbTable: prob, RevCompEnergyTable: energyRevComp, RevCompProbTable: probRevComp}
}

// Temporario: imprime os valores de energia e probabilidade dos cromossomos em csv
// func print(s []float32) {
// 	fmt.Println("seqname\tstart\tend\tscore")
// 	l := 8
// 	for i := 0; i < len(s); i++ {
// 		if s[i] > 0.5 {
// 			fmt.Printf("chrY\t%d\t%d\t%.2f\n", i, i+l-1, s[i])
// 		}
// 	}
// }

func (cr *ChrRec) Output() {
	fmt.Println("seqname\tstart\tend\tscore\tstrand")
	lrec := cr.RecLength
	lscore := len(cr.SenseProb)
	for i := 0; i < lscore; i++ {
		if cr.SenseProb[i] > 0.5 {
			fmt.Printf("%s\t%d\t%d\t%.2f\t+\n", cr.ChrName, i+1, i+lrec, cr.SenseProb[i])
		}
		if cr.AntiSenseProb[i] > 0.5 {
			fmt.Printf("%s\t%d\t%d\t%.2f\t-\n", cr.ChrName, i+1, i+lrec, cr.AntiSenseProb[i])
		}
	}

}

func main() {
	fnchr := flag.String("chr", "", "chromossome input fasta file")
	fntable := flag.String("table", "", "receptor - energy table (IMPORTANT: descending order)")

	flag.Parse()

	if *fnchr == "" || *fntable == "" {
		fmt.Println("Please use '-h' to see usage")
		return
	}

	chName, chSeq := readChr(*fnchr)
	rec := readTable("RECEPTOR NAME", *fntable)
	chrRec := rec.Scan(chName, chSeq)
	chrRec.Output()

	// print(chrRecA.SenseProb)

	// ch := readChr("/home/jgcarvalho/gocode/src/github.com/jgcarvalho/TFBS-scan/Hs.chrY.fasta")

	// recA := readTable("P10589", "/home/jgcarvalho/gocode/src/github.com/jgcarvalho/TFBS-scan/P10589.csv")
	// recB := readTable("P19793", "/home/jgcarvalho/gocode/src/github.com/jgcarvalho/TFBS-scan/P19793.csv")
	// chrRecA := recA.Scan("chrY", ch)
	// chrRecB := recB.Scan("chrY", ch)
	// fmt.Println("Done")
	// fmt.Println(scoreEnergy[100000:100100])

	// fmt.Println(chrRecA.SenseProb[100000:100100])
	// fmt.Println(chrRecB.SenseProb[100000:100100])

	// fmt.Println(chrData.SenseEnergy[100000:100100])
	// fmt.Println(chrData.SenseProb[100000:100100])
	// fmt.Println(chrData.AntiSenseEnergy[100000:100100])
	// fmt.Println(chrData.AntiSenseProb[100000:100100])
	// fmt.Println(chrRecA.SenseProb)

	// print(chrRecA.SenseProb)

}
