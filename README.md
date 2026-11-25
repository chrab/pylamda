# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/chrab/pylamda/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                        |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|---------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| src/pylamda/\_\_init\_\_.py |       16 |        2 |        0 |        0 |     88% |     49-50 |
| src/pylamda/\_\_main\_\_.py |        3 |        3 |        0 |        0 |      0% |       3-7 |
| src/pylamda/bib.py          |      483 |       37 |      234 |       15 |     92% |136-\>135, 150-\>153, 184-\>179, 214, 221, 375, 392-393, 399-\>exit, 422-430, 435, 438-447, 451-452, 456, 460-467, 706, 747-\>751, 768-769, 771-\>773, 787-\>783, 812-\>807, 850-\>854, 856, 858 |
| src/pylamda/catalogue.py    |      116 |        3 |       30 |        1 |     97% |111, 281-282 |
| src/pylamda/cli.py          |      388 |       20 |      116 |        5 |     95% |310-312, 472-474, 484-489, 491-\>493, 505-507, 510-511, 527-529, 537, 636-\>633 |
| src/pylamda/combine.py      |      268 |        6 |      104 |        6 |     97% |211, 298, 301, 400, 437-\>439, 684-685 |
| src/pylamda/compare.py      |      449 |       11 |      160 |       18 |     95% |152, 191, 381-\>379, 608-\>610, 611, 613, 625, 627, 628-\>638, 633-\>638, 639, 667, 689, 777, 900, 901-\>910, 1002-\>1006, 1009-\>1018 |
| src/pylamda/fix.py          |      443 |       18 |      220 |       16 |     95% |158, 161, 197-\>201, 203, 241-\>234, 328-\>325, 334-\>331, 356, 401, 438, 457, 512, 538-\>537, 605, 645, 662, 824-830 |
| src/pylamda/formula.py      |       36 |        2 |       16 |        2 |     92% |    64, 73 |
| src/pylamda/model.py        |      139 |        0 |       22 |        3 |     98% |440-\>442, 442-\>451, 452-\>451 |
| src/pylamda/reader.py       |      272 |        3 |       98 |        6 |     98% |209, 215, 325-\>327, 351-\>353, 353-\>355, 432 |
| src/pylamda/references.py   |       44 |        0 |       22 |        4 |     94% |64-\>62, 67-\>66, 88-\>87, 114-\>116 |
| src/pylamda/refsearch.py    |      372 |       76 |      162 |       21 |     76% |187-202, 217-\>208, 229-\>231, 234, 242, 285, 310-\>312, 323-\>326, 331-337, 343, 369, 375, 383, 418-419, 421-422, 430, 438-439, 443-446, 449-450, 454, 473, 479-\>exit, 504-529, 605-612 |
| src/pylamda/verify.py       |      372 |        4 |      156 |        7 |     98% |182, 249, 449-\>467, 507-\>509, 641, 644, 661-\>632 |
| src/pylamda/writer.py       |       77 |        0 |       36 |        2 |     98% |136-\>140, 138-\>140 |
| **TOTAL**                   | **3514** |  **185** | **1376** |  **106** | **93%** |           |

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