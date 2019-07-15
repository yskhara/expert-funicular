# expert-funicular
trajectory generation

軌道生成しまーす
擬スペクトル法によって直接最適化しまーす
最適化に関する数学操作はすべて，優れた最適化ツールであるIPOPTとPSOPTがやってくれまーす
デフォルトだとルジャンドル多項式によって内点を定めるルジャンドル擬スペクトル法を使用しまーす
[https://github.com/coin-or/Ipopt](IPOPT)と[PSOPT](https://github.com/PSOPT/psopt)などが必要でーす

PSOPTをインストールすると勝手にIPOPTもインストールされると思いまーす

世に出回っている最適化ツールの殆どはプロプラエタリですが，IPOPTはオープンで，誰でもが使用でき，優れたパフォーマンスを示します．
IPOPTやPSOPT，それらが依存するコードにcontributeされた方々に感謝と敬意を表します．
