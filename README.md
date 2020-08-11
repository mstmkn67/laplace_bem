# laplace_bem

2次元,３次元のラプラス方程式の境界要素法および高速多重極展開境界要素法のプログラムです。
[OCTA](http://octa.jp/)のインターフェースであるGOURMETで動きます。
GOURMETのメニューから、 Tool>Environments Setup から、ダウンロードしたディレクトリを指定して使ってください。

実行ファイルは、Ubuntu18.04.5でコンパイルしています。 
通常の境界要素法で実行する場合は、
laplace2d -I inputfile -O outputfile で実行してください。
高速多重極展開の境界要素法で実行する場合は、
laplace_fmm2d -I inputfile -O outputfile で実行してください。
3次元の場合は、laplace3d, laplace_fmm3dとなります。
たとえば、 laplace2d -I sphere2d.udf -O sphere2d_o.udf のようにです。

-------------
高速多重極展開は、２次元なら４分木、３次元なら８分木を作って、多重極で、要素のデータを束ねて処理するので、大規模な計算が出来ます。
![tree1](tree1.png)
![tree2](tree2.png)
-----------
![dispersion](dispersion.png)
![evaporation](evaporation.png)
