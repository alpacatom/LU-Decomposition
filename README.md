LU-Decompositision
====

LU分解のガウス消去法による解法をFortranとMPIで実装しました。
実行にはmpif.hが必要ですが、無くしたので見つけ次第追加します。

## lu_decomp_before_opt.f
並列化・最適化前のプログラム。

## lu_decomp_after_opt.f
並列化・最適化後のプログラム。
前進代入と後退代入のステップでは、ウェーブフロント処理を実装しています。
