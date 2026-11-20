# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/chrab/pylamda/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                        |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|---------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| src/pylamda/\_\_init\_\_.py |       16 |        2 |        0 |        0 |     88% |     40-41 |
| src/pylamda/\_\_main\_\_.py |        3 |        3 |        0 |        0 |      0% |       3-7 |
| src/pylamda/bib.py          |      484 |       37 |      234 |       15 |     92% |136-\>135, 150-\>153, 185-\>180, 215, 222, 376, 393-394, 400-\>exit, 423-431, 436, 439-448, 452-453, 457, 461-468, 707, 748-\>752, 769-770, 772-\>774, 788-\>784, 813-\>808, 851-\>855, 857, 859 |
| src/pylamda/catalogue.py    |      116 |        3 |       30 |        1 |     97% |111, 281-282 |
| src/pylamda/cli.py          |      308 |       22 |       88 |        5 |     93% |221, 227-229, 361-363, 373-378, 380-\>382, 394-396, 399-402, 420-422, 430 |
| src/pylamda/combine.py      |      110 |        4 |       40 |        4 |     95% |118, 205, 208, 310 |
| src/pylamda/compare.py      |      442 |       11 |      160 |       18 |     95% |153, 192, 373-\>371, 595-\>597, 598, 600, 612, 614, 615-\>622, 620-\>622, 623, 651, 668, 765, 887, 888-\>897, 991-\>994, 997-\>1006 |
| src/pylamda/fix.py          |      442 |       18 |      218 |       16 |     95% |154, 191-\>195, 197, 235-\>228, 287, 322-\>319, 328-\>325, 350, 395, 432, 451, 506, 532-\>531, 599, 639, 656, 818-824 |
| src/pylamda/formula.py      |       36 |        2 |       16 |        2 |     92% |    64, 73 |
| src/pylamda/model.py        |      135 |        0 |       22 |        3 |     98% |416-\>418, 418-\>427, 428-\>427 |
| src/pylamda/reader.py       |      196 |        1 |       64 |        4 |     98% |212-\>214, 237-\>239, 239-\>241, 310 |
| src/pylamda/references.py   |       44 |        0 |       22 |        4 |     94% |64-\>62, 67-\>66, 88-\>87, 114-\>116 |
| src/pylamda/refsearch.py    |      372 |       76 |      162 |       21 |     76% |187-202, 217-\>208, 229-\>231, 234, 242, 285, 310-\>312, 323-\>326, 331-337, 343, 369, 375, 383, 418-419, 421-422, 430, 438-439, 443-446, 449-450, 454, 473, 479-\>exit, 504-529, 605-612 |
| src/pylamda/verify.py       |      347 |        2 |      140 |        3 |     99% |151, 218, 428-\>446 |
| src/pylamda/writer.py       |       67 |        1 |       34 |        4 |     95% |27, 51-\>53, 113-\>117, 115-\>117 |
| **TOTAL**                   | **3147** |  **182** | **1230** |  **100** | **93%** |           |

1 file skipped due to complete coverage.


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/chrab/pylamda/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/chrab/pylamda/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/chrab/pylamda/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/chrab/pylamda/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2Fchrab%2Fpylamda%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/chrab/pylamda/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.