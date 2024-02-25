# はじめに
**安全に資する使い方をしてください．**

このリポジトリは弾性翼のねじり弾性モード抽出，ダイバージェンス速度の計算，フゴイド連成ダイバージェンス速度の計算に対応しています．

あくまで特定の現象のみを計算したものであるため，ここで計算される速度以下で飛行すれば必ず安全なわけではないことをご留意ください．

# 導入方法
まずpythonの動く環境を用意してください．
次に，必要なモジュールをインストールするために，requirements.txtを利用してください．
使い方は以下のとおりです．
```
 pip install -r requirements.txt
```
これによって必要なモジュールはすべてインストールされ，pythonプログラムが実行可能になります．

# 各プログラムの紹介

torsionalModeDetection.py

このプログラムはwing.csvからねじり剛性・区間あたりの重量・コード長を読み込んで有限要素法でねじり弾性のモードを計算します．
計算結果は-modalのつくプログラムで使用されます．

divergence.py

このプログラムはwing.csvからねじり剛性・桁位置・コード長を読み込んで有限要素法でダイバージェンス速度を計算します．
計算の結果得られる速度に対する最大固有値プロットのゼロクロス速度がダイバージェンス速度になります．

phugoid-divergence-modal.py

これは計算した弾性モードに基づいてフゴイド運動と連成した場合のダイバージェス速度を計算します．
当初は有限要素法で計算を行っていましたが，この場合桁弾性に150次元，フゴイド運動に1次元と非常に行列の成分数のバランスが悪くなり，数値計算上の不都合が生じたためモードでの計算を行いました．
divergence.pyと同様に，計算の結果得られる速度に対する最大固有値プロットのゼロクロス速度がダイバージェンス速度になります．

# wing.csvの構成
wing.csvはヘッダー行とそれに対応する値が列として入っている，pandasで使用する前提のcsvファイルです．
wing.csvを書き換える場合，まずspanにmm単位の翼根から翼端までの代表点位置を記述して，その代表点位置の示す行に，代表点位置におけるねじり剛性などを入力します．
不等間隔のデータでも補完して等間隔に修正されて使用されます．これはscipy.interpolate.interp1dの許容するデータであれば問題ないです．
