$nmin = 16
$nmax = 19

for ($i=$nmin; $i -le $nmax; $i++) {
    $fname = 'aln_seqs_run' + $i + '_exp' + $i + '.fa'
    $wordname = 'words_' + $i '.txt'
    
    python D:\Google Drive\To File-Research Docs\Screening\nextgenAnalysis\nextgen4B\ngsMultiExtract.py $fname $wordname 30 45 61 76 91
}
