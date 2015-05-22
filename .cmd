set terminal postscript eps enhanced dashed dashlength 1.0 linewidth 1.0 size 3.5,3 font 'Calibri,14' fontfile 'calibri.pfb' fontfile 'GillSansMT.pfb' fontfile 'GillSansItalic.pfb'
set out '.eps'
set title 'Allele Frequency Spectrum' font 'GillSansMT,18'
set grid x y mx my
set logscale xy
set key below box samplen 1 width -2
set xtics nomirror (1,2,3,4,5,10,20,50,100,200,500,1000,2000,5000)out
set ytics nomirror (1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5) out
set xlabel 'Non-reference allele count (AC)'
set ylabel 'Fraction of variants'
plot '.dat' u 1:2 lc rgbcolor 'black' lt 1 lw 3 with lines title 'Neutral Expectation at Constant Population Size'
