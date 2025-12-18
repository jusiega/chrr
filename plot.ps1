#ok robimy taki prosty system że ja podaję nazwe folderu z katalogu głównego gdzie są dane tzw liczbowe a on mi robi rysunki w główny/fig/nazwa/
$xxx=$args[0]
mkdir "./fig/$xxx"
cp $xxx/E0.txt fig/$xxx/E0.txt
cp $xxx/info.txt fig/$xxx/info.txt
cp $xxx/split.txt fig/$xxx/split.txt
cp $xxx/fxpl.txt fig/$xxx/fxpl.txt
cp $xxx/villainousinput.hxx fig/$xxx/villainousinput.hxx
gnuplot -e "cat='$xxx'" plotmaps3.gnu 
gnuplot -e "cat='$xxx'" plot.gnu &
gnuplot -e "cat='$xxx'" plot2.gnu &
gnuplot -e "cat='$xxx'" plotmaps1.gnu &
gnuplot -e "cat='$xxx'" plotmaps1.5.gnu &
gnuplot -e "cat='$xxx'" plotmaps2.gnu &
gnuplot -e "cat='$xxx'" plotmaps2.5.gnu &
gnuplot -e "cat='$xxx'" plotmaps4.gnu &
gnuplot -e "cat='$xxx'" plottrans.gnu &
Get-Job | Wait-Job | Receive-Job 
echo "FIN ALL"
pdflatex -interaction=batchmode  -output-directory="fig/$xxx" "\def\folder{fig/$xxx}\input{plot.tex}"
