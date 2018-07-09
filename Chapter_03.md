Chapter 03. The R Programming Language
================
A Solomon Kurz
2018-07-09

3.2. A simple example of R in action
------------------------------------

Basic arithmetic is straightforward in R.

``` r
2 + 3
```

    ## [1] 5

Algebra is simple, too.

``` r
x <- 2

x + x
```

    ## [1] 4

Figure 3.1

``` r
library(tidyverse)

d <-
  tibble(x = seq(from = -2, to = 2, by = .1)) %>%
  mutate(y = x^2) 
  
ggplot(data = d,
       aes(x = x, y = y)) +
  geom_line(color = "skyblue") +
  theme(panel.grid = element_blank())
```

![](Chapter_03_files/figure-markdown_github/unnamed-chunk-3-1.png)

If you’re new to the tidyverse and/or making figures with ggplot2, it’s worthwhile to walk that code out. With the first line, `library(tidyverse)`, we opened up the [core packages within the tidyverse](https://www.tidyverse.org/packages/), which are: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, and forcats.

With the next block

``` r
d <-
  tibble(x = seq(from = -2, to = 2, by = .1)) %>%
  mutate(y = x^2) 
```

we made our tibble. In R, data frames are one of the primary types of data objects (see subsection 3.4.4., below). We'll make extensive use of data frames in this project. Tibbles are a particular type of data frame, which you might learn more about [here](http://r4ds.had.co.nz/tibbles.html). With those first two lines, we determined what the name of our tibble would be, `d`, and made the first column, `x`.

Note the `%>%` operator at the end of the second line. In pose, we call that the pipe. As explained in [chapter 5 of R4DS](http://r4ds.had.co.nz/transform.html#combining-multiple-operations-with-the-pipe), "a good way to pronounce `%>%` when reading code is 'then.'" So in words, the those first two lines indicate "Make an object, `d`, which is a tibble with a variable, `x`, defined by the `seq()` function, *then*..."

In the portion after *then* (i.e., the `%>%`), we changed `d`. The `mutate()` function let us add another variable, `y`, which is a function of our first variable, `x`.

With the next 4 lines of code, we made our plot. When plotting with ggplot2, the first line is always with the `ggplot()` function. This is where you typically tell ggplot2 what data object you’re using—which must be a data frame or tibble—and what variables you want on your axes. The interesting thing about ggplot2 is that the code is modular. So if we only coded the `ggplot()` portion, we’d get:

``` r
ggplot(data = d,
       aes(x = x, y = y))
```

![](Chapter_03_files/figure-markdown_github/unnamed-chunk-5-1.png)

Although ggplot2 knows which variables to put on which axes, it has no idea how we'd like to express the data. The result is an empty coordinate system. The next line of code is the main event. With `geom_line()` we told ggplot2 to connect the data points with a line. With the `color` argument, we made that line `skyblue`. \[[Here’s a great list](http://sape.inf.usi.ch/quick-reference/ggplot2/colour) of the named colors available in ggplot2.\] Also, notice the `+` operator at the end of the `ggplot()` function. With ggplot2, you add functions by placing the `+` operator on the right of the end of one function, which will then append the next function.

``` r
ggplot(data = d,
       aes(x = x, y = y)) +
  geom_line(color = "skyblue")
```

![](Chapter_03_files/figure-markdown_github/unnamed-chunk-6-1.png)

Personally, I’m not a fan of gridlines. They occasionally have their place and I do use them from time to time. But on the while, I prefer to omit them from my plots. The final `theme()` function allowed me to do so.

``` r
ggplot(data = d,
       aes(x = x, y = y)) +
  geom_line(color = "skyblue") +
  theme(panel.grid = element_blank())
```

![](Chapter_03_files/figure-markdown_github/unnamed-chunk-7-1.png)

[Chapter 3 of R4DS](http://r4ds.had.co.nz/data-visualisation.html) is a great introduction to plotting with ggplot2. If you want to dive deeper, see the [references at the bottom of this page](https://ggplot2.tidyverse.org).

3.3. Basic commands and operators in R
--------------------------------------

In addition to the resource link Kruschke provided in the text, Grolemund and Wickham's [R4DS](http://r4ds.had.co.nz) is an excellent general introduction to the kinds of R functions you'll want to succeed with your data analysis. Other than that, I’ve learned the most when I had a specific data problem to solve and then sought out the specific code/techniques required to solve it. If already have your own data or can get your hands on some sexy data, learn these techniques by playing around with them. This isn’t the time to worry about rigor, preregistration, or all of that. This is time to play.

### 3.3.1. Getting help in R.

As with `plot()` you can learn more about the `ggplot()` function with `?`.

``` r
?ggplot
```

`help.start()` can be nice, too.

``` r
help.start()
```

`??geom_line()` can help us learn more about the `geom_line()` function.

``` r
??geom_line()
```

### 3.3.2. Arithmetic and logical operators.

With arithmetic, the order of operations is: power first, then multiplication, then addition.

``` r
1 + 2 * 3^2
```

    ## [1] 19

Whth parentheses, you can force addition before multiplication.

``` r
(1 + 2) * 3^2
```

    ## [1] 27

Operations inside parentheses get done before power operations.

``` r
(1 + 2 * 3)^2
```

    ## [1] 49

One can nest parentheses.

``` r
((1 + 2) * 3)^2
```

    ## [1] 81

``` r
?Syntax
```

We can use R to perform a variety of logical tests, such as negation.

``` r
!TRUE
```

    ## [1] FALSE

We can do conjunction.

``` r
TRUE & FALSE
```

    ## [1] FALSE

And we can do disjunction.

``` r
TRUE | FALSE
```

    ## [1] TRUE

Conjunction has precedence over disjunction.

``` r
TRUE | TRUE & FALSE
```

    ## [1] TRUE

However, with parentheses we can force disjunction first.

``` r
(TRUE | TRUE) & FALSE
```

    ## [1] FALSE

### 3.3.3. Assignment, relational operators, and tests of equality.

In contrast to Kruschke's preference, I will use the [arrow operator, `<-`, to assign](http://style.tidyverse.org/syntax.html#assignment) values to named variables.

``` r
x = 1

x <- 1
```

Yep, this ain't normal math.

``` r
(x = 1)
```

    ## [1] 1

``` r
(x = x + 1)
```

    ## [1] 2

Here we use `==` to test for equality.

``` r
(x = 2)
```

    ## [1] 2

``` r
x == 2
```

    ## [1] TRUE

Using `!=`, we can check whether the value of `x` is NOT equal to 3.

``` r
x != 3
```

    ## [1] TRUE

We can use `<` to check whether the value of `x` is less than 3.

``` r
x < 3
```

    ## [1] TRUE

Similarly, we can use `>` to check whether the value of `x` is greater than 3.

``` r
x > 3
```

    ## [1] FALSE

This normal use of the `<-` operator

``` r
x <- 3
```

is not the same as

``` r
x < - 3
```

    ## [1] FALSE

The limited precision of a computer's memory can lead to odd results.

``` r
x <- 0.5 - 0.3

y <- 0.3 - 0.1
```

Although mathematically `TRUE`, this is `FALSE` for limited precision.

``` r
x == y
```

    ## [1] FALSE

However, they are equal up to the precision of a computer.

``` r
all.equal(x, y)
```

    ## [1] TRUE

3.4. Variable types
-------------------

### 3.4.1. Vector.

#### 3.4.1.1. The combine function.

The combine function is `c()`.

``` r
c(2.718, 3.14, 1.414)
```

    ## [1] 2.718 3.140 1.414

``` r
x <- c(2.718, 3.14, 1.414)
```

You'll note the equivalence.

``` r
x == c(2.718, 3.14, 1.414)
```

    ## [1] TRUE TRUE TRUE

This leads to the next subsection.

#### 3.4.1.2. Component-by-component vector operations.

We can multiply two vectors, component by component.

``` r
c(1, 2, 3) * c(7, 6, 5)
```

    ## [1]  7 12 15

If you have a sole number, a scaler, you can multiply an entire vector by it like:

``` r
2 * c(1, 2, 3)
```

    ## [1] 2 4 6

Which is a more compact way to perform:

``` r
c(2, 2, 2) * c(1, 2, 3)
```

    ## [1] 2 4 6

The same sensibilities hold for other operations, such as addition.

``` r
2 + c(1, 2, 3)
```

    ## [1] 3 4 5

#### 3.4.1.3. The colon operator and sequence function

The colon operator has precedence over addition.

``` r
2 + 3:6
```

    ## [1] 5 6 7 8

Parentheses override default precedence.

``` r
(2 + 3):6 
```

    ## [1] 5 6

The power operator has precedence over the colon operator.

``` r
1:3^2
```

    ## [1] 1 2 3 4 5 6 7 8 9

And parentheses override default precedence.

``` r
(1:3)^2
```

    ## [1] 1 4 9

The `seq()` function is quite handy. If you don't specify the length of the output, it will figure that out the logical consequence of the other arguments.

``` r
seq(from = 0, to = 3, by = 0.5)
```

    ## [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0

This sequence won't exceed `to = 3`.

``` r
seq(from = 0, to = 3, by = 0.5001) 
```

    ## [1] 0.0000 0.5001 1.0002 1.5003 2.0004 2.5005

In each of the following examples, we'll omit one of the core `seq()` arguments: `from`, `to`, `by`, and `length.out`. Here we do not define the end point.

``` r
seq(from = 0, by = 0.5, length.out = 7)
```

    ## [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0

This time we fail to define the increment.

``` r
seq(from = 0, to = 3, length.out = 7)
```

    ## [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0

And this time we omit a starting point.

``` r
seq(to = 3, by = 0.5, length.out = 7)
```

    ## [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0

#### 3.4.1.4. The replicate function.

We'll define our pre-replication vector with the `<-` operator.

``` r
abc <- c("A", "B", "C")
```

With `times`, we repeat the vector as a unit.

``` r
rep(abc, times = 2)
```

    ## [1] "A" "B" "C" "A" "B" "C"

But if we mix `times` with `c()`, we can repeat individual components of `abc` differently.

``` r
rep(abc, times = c(4, 2, 1))
```

    ## [1] "A" "A" "A" "A" "B" "B" "C"

With the `each` argument, we repeat the individual components of `abc` one at a time.

``` r
rep(abc, each = 2)
```

    ## [1] "A" "A" "B" "B" "C" "C"

And you can even combine `each` and `length`, repeating each element until the `length` requirement has been fulfilled.

``` r
rep(abc, each = 2, length = 10)
```

    ##  [1] "A" "A" "B" "B" "C" "C" "A" "A" "B" "B"

You can also combine `each` and `times`.

``` r
rep(abc, each = 2, times = 3)
```

    ##  [1] "A" "A" "B" "B" "C" "C" "A" "A" "B" "B" "C" "C" "A" "A" "B" "B" "C"
    ## [18] "C"

I tend to do things like the above as two separate steps. One way to do so is by nesting one `rep()` function within another.

``` r
rep(rep(abc, each = 2),
    times = 3)
```

    ##  [1] "A" "A" "B" "B" "C" "C" "A" "A" "B" "B" "C" "C" "A" "A" "B" "B" "C"
    ## [18] "C"

As Kruschke points out, this can look confusing.

``` r
rep(abc, each = 2, times = c(1, 2, 3, 1, 2, 3))
```

    ##  [1] "A" "A" "A" "B" "B" "B" "B" "C" "C" "C" "C" "C"

But breaking the results up into two steps might be easier to understand,

``` r
rep(rep(abc, each = 2), 
    times = c(1, 2, 3, 1, 2, 3))
```

    ##  [1] "A" "A" "A" "B" "B" "B" "B" "C" "C" "C" "C" "C"

And especially earlier in my R career, it helped quite a bit to break operation sequences like this up by saving and assessing the intermediary steps.

``` r
step_1 <- rep(abc, each = 2)
step_1
```

    ## [1] "A" "A" "B" "B" "C" "C"

``` r
rep(step_1, times = c(1, 2, 3, 1, 2, 3))
```

    ##  [1] "A" "A" "A" "B" "B" "B" "B" "C" "C" "C" "C" "C"

#### 3.4.1.5. Getting at elements of a vector.

Behold our exemplar vector, `x`.

``` r
x <- c(2.718, 3.14, 1.414, 47405)
```

The straightforward way to extract the second and fourth elements is

``` r
x[c(2, 4)]
```

    ## [1]     3.14 47405.00

Or you might use reverse logic and omit the first and third elements.

``` r
x[c(-1, -3 )]
```

    ## [1]     3.14 47405.00

It's handy to know that `T` is a stand in for `TRUE` and `F` is a stand in for `FALSE`.

``` r
x[c(F, T, F, T)]
```

    ## [1]     3.14 47405.00

The `names()` function makes it easy to name the components of a vector.

``` r
names(x) <- c("e", "pi", "sqrt2", "zipcode")

x
```

    ##         e        pi     sqrt2   zipcode 
    ##     2.718     3.140     1.414 47405.000

Now we can call the components with their names.

``` r
x[c("pi", "zipcode")]
```

    ##       pi  zipcode 
    ##     3.14 47405.00

Here's Kruschke's review:

``` r
# define a vector
x <- c(2.718, 3.14, 1.414, 47405)

# name the components
names(x) <- c("e", "pi", "sqrt2", "zipcode")

# you can indicate which elements you'd like to include
x[c(2, 4)]
```

    ##       pi  zipcode 
    ##     3.14 47405.00

``` r
# you can decide which to exclude

x[c(-1, -3)]
```

    ##       pi  zipcode 
    ##     3.14 47405.00

``` r
# or you can use logical tests
x[c(F, T, F, T)]
```

    ##       pi  zipcode 
    ##     3.14 47405.00

``` r
# and you can use the names themselves
x[c("pi", "zipcode")]
```

    ##       pi  zipcode 
    ##     3.14 47405.00

### 3.4.2. Factor.

Here are our five-person SES status data.

``` r
x <- c("high", "medium", "low", "high", "medium")
x
```

    ## [1] "high"   "medium" "low"    "high"   "medium"

The `factor()` function turns them into a factor, which will return the levels when called.

``` r
xf <- factor(x)
xf
```

    ## [1] high   medium low    high   medium
    ## Levels: high low medium

Here are the factor levels as numerals.

``` r
as.numeric(xf)
```

    ## [1] 1 3 2 1 3

With the `levels` and `ordered` arguments, we can order the factor elements.

``` r
xfo <- factor(x, levels = c("low", "medium", "high"), ordered = T)
xfo
```

    ## [1] high   medium low    high   medium
    ## Levels: low < medium < high

Now "high" is a larger integer.

``` r
as.numeric(xfo)
```

    ## [1] 3 2 1 3 2

We've already specified `xf`.

``` r
xf
```

    ## [1] high   medium low    high   medium
    ## Levels: high low medium

And we know how it's been coded numerically.

``` r
as.numeric(xf)
```

    ## [1] 1 3 2 1 3

We can have `levels` and `labels`.

``` r
xfol <- factor(x, 
               levels = c("low", "medium", "high"), ordered = T,
               labels = c("Bottom SES", "Middle SES", "Top SES"))

xfol
```

    ## [1] Top SES    Middle SES Bottom SES Top SES    Middle SES
    ## Levels: Bottom SES < Middle SES < Top SES

### 3.4.3. Matrix and array.

``` r
matrix(1:6, ncol = 3)
```

    ##      [,1] [,2] [,3]
    ## [1,]    1    3    5
    ## [2,]    2    4    6

We can get the same thing using `nrow`.

``` r
matrix(1:6, nrow = 2)
```

    ##      [,1] [,2] [,3]
    ## [1,]    1    3    5
    ## [2,]    2    4    6

Note how the numbers got ordered by rows within each column? We can specify them to be ordered across columns, first.

``` r
matrix(1:6, nrow = 2, byrow = T)
```

    ##      [,1] [,2] [,3]
    ## [1,]    1    2    3
    ## [2,]    4    5    6

We can name the dimensions. I'm not completely consistent, but I've been moving in the direction of following [The Tidyverse Style Guide](http://style.tidyverse.org) for naming my R objects and their elements. From [the guide](http://style.tidyverse.org/syntax.html), we read

> Variable and function names should use only lowercase letters, numbers, and `_`. Use underscores (`_`) (so called snake case) to separate words within a name.

By those sensibilities, we'll name our rows and columns as

``` r
matrix(1:6, 
       nrow = 2,
       dimnames = list(TheRowDimName = c("row_1_name", "row_2_name"),
                      TheColDimName  = c("col_1_name", "col_2_name", "col_3_name")))
```

    ##              TheColDimName
    ## TheRowDimName col_1_name col_2_name col_3_name
    ##    row_1_name          1          3          5
    ##    row_2_name          2          4          6

You've also probably noticed that I "[always put a space after a comma, never before, just like in regular English](http://style.tidyverse.org/syntax.html#spacing)," as well as "put a space before and after `=` when naming arguments in function calls." IMO, this makes code easier to read.

We'll name our matrix `x`.

``` r
x <-
  matrix(1:6, 
       nrow = 2,
       dimnames = list(TheRowDimName = c("row_1_name", "row_2_name"),
                       TheColDimName  = c("col_1_name", "col_2_name", "col_3_name")))
```

Since there are 2 dimensions, we'll subset with two dimensions. Numerical indices work.

``` r
x[2, 3]
```

    ## [1] 6

Row and column names work, too. Just make sure to use quotation marks, `""`, for those.

``` r
x["row_2_name", "col_3_name"]
```

    ## [1] 6

Here we specify the range of columns to include.

``` r
x[2, 1:3]
```

    ## col_1_name col_2_name col_3_name 
    ##          2          4          6

Leaving that argument blank returns them all.

``` r
x[2, ]
```

    ## col_1_name col_2_name col_3_name 
    ##          2          4          6

And leaving the row index blank returns all row values within the specified column(s).

``` r
x[, 3]
```

    ## row_1_name row_2_name 
    ##          5          6

Mind your commas! This produces the second row, returned as a vector.

``` r
x[2, ]
```

    ## col_1_name col_2_name col_3_name 
    ##          2          4          6

This returns both rows of the 2<sup>nd</sup> column.

``` r
x[, 2]
```

    ## row_1_name row_2_name 
    ##          3          4

Leaving out the comma will return the numbered element.

``` r
x[2]
```

    ## [1] 2

It'll be important in your brms career to have a sense of 3-dimensional arrays. Several brms convenience functions often return them (e.g., `ranef()` in multilevel models).

``` r
a <- array(1:24, dim = c(3, 4, 2), # 3 rows, 4 columns, 2 layers
           dimnames = list(RowDimName = c("r1", "r2", "r3"),
                           ColDimName = c("c1", "c2", "c3", "c4"),
                           LayDimName = c("l1", "l2")))

a
```

    ## , , LayDimName = l1
    ## 
    ##           ColDimName
    ## RowDimName c1 c2 c3 c4
    ##         r1  1  4  7 10
    ##         r2  2  5  8 11
    ##         r3  3  6  9 12
    ## 
    ## , , LayDimName = l2
    ## 
    ##           ColDimName
    ## RowDimName c1 c2 c3 c4
    ##         r1 13 16 19 22
    ##         r2 14 17 20 23
    ##         r3 15 18 21 24

Since these have 3 dimensions, you have to use 3-dimensional indexing. As with 2-dimensional objects, leaving the indices for a dimension blank will return all elements within that dimension. For example, this code returns all columns of `r3` and `l2`, as a vector.

``` r
a["r3", , "l2"]
```

    ## c1 c2 c3 c4 
    ## 15 18 21 24

And this code returns all layers of `r3` and `c4`, as a vector.

``` r
a["r3", "c4", ]
```

    ## l1 l2 
    ## 12 24

### 3.4.4. List and data frame.

Here's `my_list`.

``` r
my_list <- 
  list("a" = 1:3, 
       "b" = matrix(1:6, nrow = 2), 
       "c" = "Hello, world.")

my_list
```

    ## $a
    ## [1] 1 2 3
    ## 
    ## $b
    ##      [,1] [,2] [,3]
    ## [1,]    1    3    5
    ## [2,]    2    4    6
    ## 
    ## $c
    ## [1] "Hello, world."

To return the contents of the `a` portion of `my_list`, just do:

``` r
my_list$a
```

    ## [1] 1 2 3

We can index further within `a`.

``` r
my_list$a[2]
```

    ## [1] 2

To return the contents of the first item in our list with the double bracket, `[[]]`, do:

``` r
my_list[[1]]
```

    ## [1] 1 2 3

You can index further to return only the second element of the first list item.

``` r
my_list[[1]][2]
```

    ## [1] 2

But double brackets, `[][]`, are no good, here.

``` r
my_list[1][2]
```

    ## $<NA>
    ## NULL

To learn more, Jenny Bryan has a [great talk](https://www.youtube.com/watch?v=4MfUCX_KpdE&t=615s&frags=pl%2Cwn) discussing the role of lists within data wrangling. But here's a data frame.

``` r
d <- 
  data.frame(integers = 1:3, 
             number_names = c("one", "two", "three"))

d
```

    ##   integers number_names
    ## 1        1          one
    ## 2        2          two
    ## 3        3        three

With data frames, we can continue indexing with the `$` operator.

``` r
d$number_names
```

    ## [1] one   two   three
    ## Levels: one three two

We can also use the double bracket.

``` r
d[[2]]
```

    ## [1] one   two   three
    ## Levels: one three two

Notice how the single bracket with no comma indexes columns rather than rows.

``` r
d[2]
```

    ##   number_names
    ## 1          one
    ## 2          two
    ## 3        three

But adding the comma returns the factor-level information when indexing columns.

``` r
d[, 2]
```

    ## [1] one   two   three
    ## Levels: one three two

It works a touch differently when indexing by row.

``` r
d[2, ]
```

    ##   integers number_names
    ## 2        2          two

Let's try with a tibble, instead.

``` r
t <-
  tibble(integers = 1:3,
         number_names = c("one", "two", "three"))

t
```

    ## # A tibble: 3 x 2
    ##   integers number_names
    ##      <int> <chr>       
    ## 1        1 one         
    ## 2        2 two         
    ## 3        3 three

One difference is that tibbles default to assigning text columns as character strings rather than factors. Another difference occurs when printing large data frames versus large tibbles. Tibbles yield more compact glimpses. For more, check out [R4DS Chapter 10](http://r4ds.had.co.nz/tibbles.html).

It’s also worthwhile pointing out that within the tidyverse, you can pull out a specific column with the `select()` function. Here we select `number_names`.

``` r
t %>% 
  select(number_names)
```

    ## # A tibble: 3 x 1
    ##   number_names
    ##   <chr>       
    ## 1 one         
    ## 2 two         
    ## 3 three

Go [here](http://r4ds.had.co.nz/transform.html#select) learn more about `select()`.

3.5. Loading and saving data
----------------------------

### 3.5.1. The ~~read.csv~~ `read_csv()` and ~~read.table~~ `read_table()` functions.

Although `read.csv()` is the default CSV reader in R, the [`read_csv()` function](https://readr.tidyverse.org/reference/read_delim.html) from the [readr package](https://readr.tidyverse.org) (i.e., one of the core tidyverse packages) is a new alternative. In [comparison to base R](http://r4ds.had.co.nz/data-import.html#compared-to-base-r)'s `read.csv()`, `readr::read_csv()` is faster and returns tibbles (as opposed to data frames with `read.csv()`). The same general points hold for base R's `read.table()` versus `readr::read_table()`.

Using Kruschke's `HGN.csv` example, we'd load the CSV with `read_csv()` like this:

``` r
hgn <- read_csv("data.R/HGN.csv")
```

Note again that `read_csv()` defaults to returning columns with character information as characters, not factors.

``` r
hgn$Hair
```

    ## [1] "black" "brown" "blond" "black" "black" "red"   "brown"

See? As a character variable, `Hair` no longer has factor level information. But if you knew you wanted to treat `Hair` as a factor, you could easily convert it with `mutate()`.

``` r
hgn <-
  hgn %>% 
  mutate(Hair = factor(Hair))
         
hgn$Hair
```

    ## [1] black brown blond black black red   brown
    ## Levels: black blond brown red

And here's a tidyverse way to reorder the levels for the `Hair` factor.

``` r
hgn <-
  hgn %>% 
  mutate(Hair = factor(Hair, levels = c("red", "blond", "brown", "black")))
         
hgn$Hair
```

    ## [1] black brown blond black black red   brown
    ## Levels: red blond brown black

``` r
as.numeric(hgn$Hair)
```

    ## [1] 4 3 2 4 4 1 3

Since we imported `hgn` with `read_csv()`, the `Name` column is already a character vector, which we can verify with the `str()` function.

``` r
hgn$Name %>% str()
```

    ##  chr [1:7] "Alex" "Betty" "Carla" "Diane" "Edward" "Frank" "Gabrielle"

Note how using `as.vector()` did nothing in our case. `Name` was already a character vector.

``` r
hgn$Name %>% 
  as.vector() %>% 
  str()
```

    ##  chr [1:7] "Alex" "Betty" "Carla" "Diane" "Edward" "Frank" "Gabrielle"

The `Group` column was imported as composed of integers.

``` r
hgn$Group %>% str()
```

    ##  int [1:7] 1 1 1 2 2 2 2

Switching `Group` to a factor is easy enough.

``` r
hgn <-
  hgn %>% 
  mutate(Group = factor(Group))

hgn$Group
```

    ## [1] 1 1 1 2 2 2 2
    ## Levels: 1 2

### 3.5.2. Saving data from R.

Yeah you guessed, readr has a `write_csv()` function, too. The arguments are as follows: `write_csv(x, path, na = "NA", append = FALSE, col_names = !append)`. Saving `hgn` in your working directory is as easy as:

``` r
write_csv(hgn, "hgn.csv")
```

You could also use `save()`.

``` r
save(hgn, file = "hgn.Rdata" )
```

Once we start fitting Bayesian models, this method will be an important way to save the results of those models.

The `load()` function is simple.

``` r
load("hgn.Rdata" )
```

The `ls()` function works very much the same way as the more verbosely-named `objects()` function.

``` r
ls()
```

    ##  [1] "a"       "abc"     "d"       "hgn"     "my_list" "step_1"  "t"      
    ##  [8] "x"       "xf"      "xfo"     "xfol"    "y"

3.6. Some utility functions
---------------------------

``` r
# This is a more compact way to replicate 100 1’s, 200 2’s, and 300 3’s
x <- rep(1:3, times = c(100, 200, 300))

summary(x)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   1.000   2.000   2.500   2.333   3.000   3.000

We can use the pipe to convert and then summarize `x`.

``` r
x %>% 
  factor() %>% 
  summary()
```

    ##   1   2   3 
    ## 100 200 300

`head()` and `tail()` are quite useful.

``` r
head(x)
```

    ## [1] 1 1 1 1 1 1

``` r
tail(x)
```

    ## [1] 3 3 3 3 3 3

Within the tidyverse, the `slice()` function serves a similar role. In order to use `slice()`, we'll want to convert `x`, which is just a vector of integers, into a data frame. Then we'll use `slice()` to return a subset of the rows.

``` r
x <-
  x %>%
  as_tibble() 

x %>% 
  slice(1:6)
```

    ## # A tibble: 6 x 1
    ##   value
    ##   <int>
    ## 1     1
    ## 2     1
    ## 3     1
    ## 4     1
    ## 5     1
    ## 6     1

So that was analogous to what we accomplished with `head()`. Here's the analogue to `tail()`.

``` r
x %>%
  slice(595:600)
```

    ## # A tibble: 6 x 1
    ##   value
    ##   <int>
    ## 1     3
    ## 2     3
    ## 3     3
    ## 4     3
    ## 5     3
    ## 6     3

The downside of that code was we had to do the math to determine that 600 - 6 = 595 in order to get the last six rows, as returned by `tail()`. A more general approach is to use `n()`, which will return the total number of rows in the tibble.

``` r
x %>%
  slice((n() - 6):n())
```

    ## # A tibble: 7 x 1
    ##   value
    ##   <int>
    ## 1     3
    ## 2     3
    ## 3     3
    ## 4     3
    ## 5     3
    ## 6     3
    ## 7     3

To unpack `(n() - 6):n()`, because `n()` = 600, `(n() - 6)` = 600 - 6 = 595. Therefore `(n() - 6):n()` was equivalent to having coded `595:600`. Instead of having to do the math ourselves, `n()` did it for us. It’s often easier to just go with `head()` or `tail()`. But the advantage of this more general approach is that it allows one take more complicated slices of the data, such as returning the first three and last three rows.

``` r
x %>%
  slice(c(1:3, (n() - 3):n()))
```

    ## # A tibble: 7 x 1
    ##   value
    ##   <int>
    ## 1     1
    ## 2     1
    ## 3     1
    ## 4     3
    ## 5     3
    ## 6     3
    ## 7     3

We've already used the handy `str()` function a bit. It's also nice to know that `tidyverse::glimpse()` performs a similar function.

``` r
x %>% str()
```

    ## Classes 'tbl_df', 'tbl' and 'data.frame':    600 obs. of  1 variable:
    ##  $ value: int  1 1 1 1 1 1 1 1 1 1 ...

``` r
x %>% glimpse()
```

    ## Observations: 600
    ## Variables: 1
    ## $ value <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,...

Within the tidyverse, we'd use `group_by()` and then `summarize()` as alternatives to the `aggregate()` function. With `group_by()` we group the observations first by `Hair` and then by `Gender` within `Hair`. After that, we summarize the groups by taking the `median()` values of their `Number`.

``` r
hgn %>% 
  group_by(Hair, Gender) %>% 
  summarize(median = median(Number))
```

    ## # A tibble: 5 x 3
    ## # Groups:   Hair [?]
    ##   Hair  Gender median
    ##   <fct> <chr>   <dbl>
    ## 1 red   M         7  
    ## 2 blond F         3  
    ## 3 brown F         7  
    ## 4 black F         7  
    ## 5 black M         1.5

One of the nice things about this workflow is that the code reads somewhat like how we’d explain what we were doing. We, in effect, told R to *Take `hgn`, then group the data by `Hair` and `Gender` within `Hair`, and then `summarize()` those groups by their `median()` `Number` values.* There’s also the nice quality that we don’t have to continually tell R where the data are coming from the way the `aggregate()` function required Kruschke to prefix each of his variables with `HGNdf$`. We also didn't have to explicitly rename the output columns the way Kruschke had to.

I'm not aware that our `group_by() %>% summarize()` workflow has a formula format the way `aggregate()` does.

To count how many levels we had in a grouping factor, we'd use the `n()` function in `summarize()`.

``` r
hgn %>% 
  group_by(Hair, Gender) %>% 
  summarize(n = n())
```

    ## # A tibble: 5 x 3
    ## # Groups:   Hair [?]
    ##   Hair  Gender     n
    ##   <fct> <chr>  <int>
    ## 1 red   M          1
    ## 2 blond F          1
    ## 3 brown F          2
    ## 4 black F          1
    ## 5 black M          2

Alternatively, we could switch out the `summary(n = n())` line with `count()`.

``` r
hgn %>% 
  group_by(Hair, Gender) %>% 
  count()
```

    ## # A tibble: 5 x 3
    ## # Groups:   Hair, Gender [5]
    ##   Hair  Gender     n
    ##   <fct> <chr>  <int>
    ## 1 red   M          1
    ## 2 blond F          1
    ## 3 brown F          2
    ## 4 black F          1
    ## 5 black M          2

We could then use `spread()` to convert that output to a format similar to Kruschke's table of counts.

``` r
hgn %>% 
  group_by(Hair, Gender) %>% 
  count() %>% 
  spread(key = Hair, value = n)
```

    ## # A tibble: 2 x 5
    ## # Groups:   Gender [2]
    ##   Gender   red blond brown black
    ##   <chr>  <int> <int> <int> <int>
    ## 1 F         NA     1     2     1
    ## 2 M          1    NA    NA     2

With this method, the `NA`s are stand-ins for 0s.

``` r
a
```

    ## , , LayDimName = l1
    ## 
    ##           ColDimName
    ## RowDimName c1 c2 c3 c4
    ##         r1  1  4  7 10
    ##         r2  2  5  8 11
    ##         r3  3  6  9 12
    ## 
    ## , , LayDimName = l2
    ## 
    ##           ColDimName
    ## RowDimName c1 c2 c3 c4
    ##         r1 13 16 19 22
    ##         r2 14 17 20 23
    ##         r3 15 18 21 24

`apply()` is part of a family of functions that offer a wide array of uses. You can learn more about the `apply()` family [here](https://www.datacamp.com/community/tutorials/r-tutorial-apply-family) or [here](http://faculty.nps.edu/sebuttre/home/R/apply.html).

``` r
apply(a, MARGIN = c(2, 3), FUN = sum)
```

    ##           LayDimName
    ## ColDimName l1 l2
    ##         c1  6 42
    ##         c2 15 51
    ##         c3 24 60
    ##         c4 33 69

Here's `a`.

``` r
a
```

    ## , , LayDimName = l1
    ## 
    ##           ColDimName
    ## RowDimName c1 c2 c3 c4
    ##         r1  1  4  7 10
    ##         r2  2  5  8 11
    ##         r3  3  6  9 12
    ## 
    ## , , LayDimName = l2
    ## 
    ##           ColDimName
    ## RowDimName c1 c2 c3 c4
    ##         r1 13 16 19 22
    ##         r2 14 17 20 23
    ##         r3 15 18 21 24

The reshape2 package [is a precursor](https://tidyr.tidyverse.org) to the tidyr package (i.e., one of the core tidyverse packages). The `reshape2::melt()` function is a quick way to transform the 3-dimensional `a` matrix into a tidy data frame.

``` r
a %>% 
  reshape2::melt()
```

    ##    RowDimName ColDimName LayDimName value
    ## 1          r1         c1         l1     1
    ## 2          r2         c1         l1     2
    ## 3          r3         c1         l1     3
    ## 4          r1         c2         l1     4
    ## 5          r2         c2         l1     5
    ## 6          r3         c2         l1     6
    ## 7          r1         c3         l1     7
    ## 8          r2         c3         l1     8
    ## 9          r3         c3         l1     9
    ## 10         r1         c4         l1    10
    ## 11         r2         c4         l1    11
    ## 12         r3         c4         l1    12
    ## 13         r1         c1         l2    13
    ## 14         r2         c1         l2    14
    ## 15         r3         c1         l2    15
    ## 16         r1         c2         l2    16
    ## 17         r2         c2         l2    17
    ## 18         r3         c2         l2    18
    ## 19         r1         c3         l2    19
    ## 20         r2         c3         l2    20
    ## 21         r3         c3         l2    21
    ## 22         r1         c4         l2    22
    ## 23         r2         c4         l2    23
    ## 24         r3         c4         l2    24

We have an alternative if you wanted to stay within the tydiverse. To my knowledge, the fastest way to make the transformation is to first use `as.tbl_cube()` and follow that up with `as_tibble()`. The [`as.tbl_cube()` function](https://dplyr.tidyverse.org/reference/as.tbl_cube.html) will convert the `a` matrix into a tbl\_cube. We will use the `met_name` argument to determine the name of the measure assessed in the data. Since the default is for `as.tbl_cube()` to name the measure name as `.`, it seemed `value` was a more descriptive choice. We'll then use the `as_tibble()` function to convert our tbl\_cube object into a tidy tibble.

``` r
a %>% 
  as.tbl_cube(met_name = "value") %>% 
  as_tibble()
```

    ## # A tibble: 24 x 4
    ##    RowDimName ColDimName LayDimName value
    ##    <chr>      <chr>      <chr>      <int>
    ##  1 r1         c1         l1             1
    ##  2 r2         c1         l1             2
    ##  3 r3         c1         l1             3
    ##  4 r1         c2         l1             4
    ##  5 r2         c2         l1             5
    ##  6 r3         c2         l1             6
    ##  7 r1         c3         l1             7
    ##  8 r2         c3         l1             8
    ##  9 r3         c3         l1             9
    ## 10 r1         c4         l1            10
    ## # ... with 14 more rows

Notice how the first three columns are returned as characters instead of factors. If you really wanted those to be factors, you could always follow up the code with `mutate_if(is.character, as.factor)`.

3.7. Programming in R
---------------------

It's worthy to note that this project was done with [R Markdown](https://rmarkdown.rstudio.com), which is an alternative to an R script. As [Grolemund and Wickham point out](http://r4ds.had.co.nz/r-markdown.html)

> R Markdown integrates a number of R packages and external tools. This means that help is, by-and-large, not available through ?. Instead, as you work through this chapter, and use R Markdown in the future, keep these resources close to hand:
>
> -   R Markdown Cheat Sheet: Help &gt; Cheatsheets &gt; R Markdown Cheat Sheet,
>
> -   R Markdown Reference Guide: Help &gt; Cheatsheets &gt; R Markdown Reference Guide.
>
> Both cheatsheets are also available at <http://rstudio.com/cheatsheets>.

I also strongly recommend checking out [R Notebooks](https://bookdown.org/yihui/rmarkdown/notebook.html), which is a kind of R Markdown document but with a few bells a whistles that make it more useful for working scientists. You can learn more about it [here](https://rstudio-pubs-static.s3.amazonaws.com/256225_63ebef4029dd40ef8e3679f6cf200a5a.html) and [here](https://www.r-bloggers.com/why-i-love-r-notebooks-2/).

### 3.7.1. Variable names in R.

Kruschke prefers to use camelBack notation for his variable and function names. Though I initially hated it, I've been moving in the direction of [snake\_case](http://style.tidyverse.org/syntax.html#object-names). If seems easier to read\_prose\_in\_snake\_case than it is to readProseInCamelBack. To each their own.

### 3.7.3. Programming a function.

Here's our simple `a_sq_plus_b` function.

``` r
a_sq_plus_b <- function(a, b = 1) {
  c <- a^2 + b
  return(c)
  }
```

If you explicitly denote your arguments, everything works fine.

``` r
a_sq_plus_b(a = 3, b = 2)
```

    ## [1] 11

Keep things explicit and you can switch up the order of the arguments.

``` r
a_sq_plus_b(b = 2, a = 3)
```

    ## [1] 11

But here's what happens when you are less explicit.

``` r
# this
a_sq_plus_b(3, 2)
```

    ## [1] 11

``` r
# is not the same as this
a_sq_plus_b(2, 3)
```

    ## [1] 7

Since we gave `b` a default value, we can be really lazy.

``` r
a_sq_plus_b(a = 2)
```

    ## [1] 5

But we can't be lazy with `a`. This

``` r
a_sq_plus_b(b = 1)
```

yielded this warning on my computer: "Error in a\_sq\_plus\_b(b = 1) : argument "a" is missing, with no default".

If we're completely lazy, `a_sq_plus_b()` presumes our sole input value is for the `a` argument and it uses the default value of `1` for `b`.

``` r
a_sq_plus_b(2)
```

    ## [1] 5

The lesson is important because it’s good practice to familiarize yourself with the defaults of the functions you use in statistics and data analysis, more generally.

### 3.7.4. Conditions and loops.

Here's our starting point for `if()` and `else()`.

``` r
if(x <= 3){      # if x is less than or equal to 3
  show("small")  # display the word "small"
  } else {       # otherwise
    show("big")  # display the word "big"
    }            # end of ’else’ clause
```

    ## Warning in if (x <= 3) {: the condition has length > 1 and only the first
    ## element will be used

    ## [1] "small"

Yep, this is no good.

``` r
if (x <= 3) {show("small")}
else {show("big")}
```

On my computer, it returned this message: "the condition has length &gt; 1 and only the first element will be used\[1\] "small" Error: unexpected 'else' in "else"".

Here we use the loop.

``` r
for (count_down in 5:1) {
  show(count_down)
  }
```

    ## [1] 5
    ## [1] 4
    ## [1] 3
    ## [1] 2
    ## [1] 1

``` r
for (note in c("do", "re", "mi")) {
  show(note)
  }
```

    ## [1] "do"
    ## [1] "re"
    ## [1] "mi"

It's also useful to understand how to use the `ifelse()` function within the context of a data frame. Recall hos `x` is a data frame.

``` r
x <- tibble(x = 1:5)

x
```

    ## # A tibble: 5 x 1
    ##       x
    ##   <int>
    ## 1     1
    ## 2     2
    ## 3     3
    ## 4     4
    ## 5     5

We can use the `mutate()` function to make a new variable, `size`, which is itself a function of the original variable, `x`. We'll use the `ifelse()` function to return "small" if `x <= 3`, but to return "big" otherwise.

``` r
x %>% 
  mutate(size = ifelse(x <= 3, "small", "big"))
```

    ## # A tibble: 5 x 2
    ##       x size 
    ##   <int> <chr>
    ## 1     1 small
    ## 2     2 small
    ## 3     3 small
    ## 4     4 big  
    ## 5     5 big

### 3.7.5. Measuring processing time.

This will be nontrivial to consider in your Bayesian career. Here's the loop.

``` r
start_time               <- proc.time()
y                        <- vector(mode = "numeric", length = 1.0E6)
for (i in 1:1.0E6) {y[i] <- log(i)}
stop_time <- proc.time()

elapsed_time_loop <- stop_time - start_time
show(elapsed_time_loop)
```

    ##    user  system elapsed 
    ##   0.101   0.008   0.116

Now we use a vector.

``` r
start_time <- proc.time()
y          <- log(1:1.0E6)
stop_time  <- proc.time()

elapsed_time_vector <- stop_time - start_time
show(elapsed_time_vector)
```

    ##    user  system elapsed 
    ##   0.021   0.007   0.029

Here we compare the two times.

``` r
elapsed_time_vector[1]/elapsed_time_loop[1]
```

    ## user.self 
    ## 0.2079208

For my computer, the vectorized approach took about 20.8% the time the loop approach did. When using R, avoid loops for vectorized approaches whenever possible.

### 3.7.6. Debugging.

This should be no surprise, by now, but in addition to Kruschke’s good advice, I also recommend checking out [R4DS](http://r4ds.had.co.nz). I reference it often.

3.8. Graphical plots: Opening and saving
----------------------------------------

For making and saving plots with ggplot2, I recommend reviewing Chapters [3](http://r4ds.had.co.nz/data-visualisation.html) and [28](http://r4ds.had.co.nz/graphics-for-communication.html) of R4DS.

References
----------

Kruschke, J. K. (2015). *Doing Bayesian data analysis, Second Edition: A tutorial with R, JAGS, and Stan.* Burlington, MA: Academic Press/Elsevier.

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.4
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2.2  forcats_0.3.0   stringr_1.3.1   dplyr_0.7.6    
    ##  [5] purrr_0.2.5     readr_1.1.1     tidyr_0.8.1     tibble_1.4.2   
    ##  [9] ggplot2_3.0.0   tidyverse_1.2.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.4 reshape2_1.4.3   haven_1.1.2      lattice_0.20-35 
    ##  [5] colorspace_1.3-2 htmltools_0.3.6  yaml_2.1.19      utf8_1.1.4      
    ##  [9] rlang_0.2.1      pillar_1.2.3     foreign_0.8-70   glue_1.2.0      
    ## [13] withr_2.1.2      modelr_0.1.2     readxl_1.1.0     bindr_0.1.1     
    ## [17] plyr_1.8.4       munsell_0.5.0    gtable_0.2.0     cellranger_1.1.0
    ## [21] rvest_0.3.2      psych_1.8.4      evaluate_0.10.1  labeling_0.3    
    ## [25] knitr_1.20       parallel_3.5.1   broom_0.4.5      Rcpp_0.12.17    
    ## [29] scales_0.5.0     backports_1.1.2  jsonlite_1.5     mnormt_1.5-5    
    ## [33] hms_0.4.2        digest_0.6.15    stringi_1.2.3    grid_3.5.1      
    ## [37] rprojroot_1.3-2  cli_1.0.0        tools_3.5.1      magrittr_1.5    
    ## [41] lazyeval_0.2.1   crayon_1.3.4     pkgconfig_2.0.1  xml2_1.2.0      
    ## [45] lubridate_1.7.4  assertthat_0.2.0 rmarkdown_1.10   httr_1.3.1      
    ## [49] rstudioapi_0.7   R6_2.2.2         nlme_3.1-137     compiler_3.5.1
