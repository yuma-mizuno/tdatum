# TDatum の公開案

2026-09-05 時点。公開対象は **SageMath で使えるライブラリと代表的な使用例**とする。

## 推奨する公開構成

独立した GitHub リポジトリを開発・配布の拠点にし、代表例をノートブックで添付する。
初回の安定版を GitHub Release として公開し、Zenodo に保存して引用用 DOI を付ける。
解説をブラウザで読めるようにする場合は、実行済みノートブックから生成した HTML を
GitHub Pages に配置する。

GitHub Pages は HTML などの静的ファイルを公開するサービスである。
SageMath の計算をサーバー側で実行する機能は含まれないため、計算は読者の SageMath
環境で行う。[GitHub Pages の公式説明](https://docs.github.com/en/pages/getting-started-with-github-pages/what-is-github-pages)

Zenodo は、接続した GitHub リポジトリの新しいリリースを取り込んで保存する。
GitHub Release の作成前に連携を有効にし、取り込み後に DOI とメタデータを確認する。
[連携の有効化](https://help.zenodo.org/docs/github/enable-repository/)・
[リリースの保存](https://help.zenodo.org/docs/github/archive-software/github-upload/)

## 公開経路の比較

| 経路 | TDatum に対する役割 | 提案 |
| --- | --- | --- |
| GitHub + Releases | ソース、変更履歴、バグ報告、配布物 | 最初の公開先 |
| GitHub Pages | 導入解説と実行済みの例を読む入口 | 初回または直後に追加 |
| Zenodo | 版を指定して引用できる保存先 | 最初の安定版と連携 |
| PyPI | Sage 環境内へのインストールを簡単にする配布先 | パッケージ名と公開 API を決めてから追加 |
| SageMath 本体への提案 | Sage の既存 API としての長期的な統合 | 外部パッケージとして運用した後に検討 |
| JOSS | 研究ソフトウェアとしての査読付き論文 | 公開開発と利用実績を蓄積した後に検討 |

現段階では外部パッケージとして公開する構成が扱いやすい。Sage では外部パッケージの
配布と本体への組み込みに別の作業があり、本体への採用を初回公開の前提にする必要はない。
[SageMath 外部パッケージ一覧](https://wiki.sagemath.org/SageMathExternalPackages)・
[パッケージ組み込みの開発者向け説明](https://doc.sagemath.org/html/en/developer/packaging.html)

JOSS の現行要件には、6 か月を超える公開開発履歴、研究での利用実績、テスト・文書化
などが含まれる。したがって初回公開と同時の投稿先としては想定せず、継続的に開発・利用
される段階で改めて適合性を判断する。
[JOSS 投稿要件](https://joss.readthedocs.io/en/latest/submitting.html)

## 今回のローカル構成

```text
tdatum/
  src/tdatum/       TDatum、MutationLoop、例の構成子
  examples/         実行済みノートブック 3 冊
  tests/            数学的性質と不正入力に対する回帰テスト
  tools/            テスト、配布物検証、ノートブック生成
  docs/             API、検証記録、出典、公開方針
  pyproject.toml    配布設定
  CITATION.cff      ソフトウェアと数学的文献の引用情報
  LICENSE          GPL 本文
```

既存の探索スクリプト、大量の探索結果、研究用の途中ノートをライブラリの導入経路に
含めない。必要になった研究例は、その入力・再現手順・結論の範囲が揃った単位で追加する。

公開用の基本導線は、README の導入手順、最小例、非自明な symmetrizer の例、RSG の例、
mutation loop から T-datum を構成する例、API 表の順とする。英語を共通の入口とし、
必要に応じて日本語の解説を併設する。

## 初回公開までの手順

1. リポジトリ名、配布名、公開 API の範囲を確定する。現状の候補はリポジトリ `tdatum`、
   配布名 `sagemath-tdatum`、Python import 名 `tdatum`。
2. 現在の `0.1.0.dev0` をレビューする。特に、旧コードからの検証条件の変更と、
   今回新たに検証していない構成子の扱いを確認する。
3. 公開用 GitHub リポジトリを作成し、ソース、ノートブック、テストを push する。
4. GitHub Actions 上でテストを実行し、SageMath 10.8 の CI 結果を確認する。
   ローカルで成功した結果と、未実行の外部 CI は区別する。
5. 安定版番号とリポジトリ URL を `pyproject.toml`、`CITATION.cff`、README に反映する。
6. Zenodo 連携を有効にし、その後にタグと GitHub Release を作成する。
7. Zenodo の DOI・著者・版・ライセンス・収録ファイルを確認し、README の引用情報に反映する。
8. 必要なら HTML 解説を GitHub Pages に配置し、次の段階で PyPI からの配布を追加する。

`CITATION.cff` は GitHub の引用表示にも使える。Zenodo 用の `.zenodo.json` を併設すると
Zenodo はそちらを優先するため、まずは `CITATION.cff` に情報を集約する。
[Zenodo の引用情報の説明](https://help.zenodo.org/docs/github/describe-software/citation-file/)

この作業ではローカルのパッケージと公開案まで作成した。公開先のアカウント操作、push、
リリース作成、DOI の発行、PyPI へのアップロードは実施していない。
