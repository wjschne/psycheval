# Convert ability (in W scores by default) to relative proficiency index

Convert ability (in W scores by default) to relative proficiency index

## Usage

``` r
rpi(
  x,
  mu = 500,
  scale = 20/log(9),
  criterion = 0.9,
  reverse = FALSE,
  interpretation = FALSE
)
```

## Arguments

- x:

  numeric vector of ability scores

- mu:

  numeric vector of ability scores of reference group

- scale:

  number vector of scaling factor. The default value (`log(9) / 20`)
  assumes that x and mu are W scores.

- criterion:

  numeric proficiency criterion (between 0 and 1, exclusive)

- reverse:

  boolean. If TRUE, the criterion refers to the proficiency of the
  person instead of the proficiency of the peer group. In other words,
  the role of the x and mu are reversed.

- interpretation:

  If TRUE, the rpi's print method will provide an interpretation of the
  relative proficiency.

## Value

numeric

## Examples

``` r
# What is the probability a person with a W score of 540 can pass
# an item that a person with a 500 W score can pass with a
# probability of .90?
rpi(x = 540, mu = 500, criterion = .9)
#> 0.9986301
# Same as above but with an interpretive statement
rpi(x = 540, mu = 500, criterion = .9, interpretation = TRUE)
#> When a same-age peer of average ability has a 1 probability of answering an item correctly, this person has a 1 probability of answering it correctly.
# When a person with a W score of 540 has a .9 probability of
# passing an item, what is the probability that a person with a W
# score of 500 will pass it?
rpi(x = 540, mu = 500, criterion = .9, reverse = TRUE, interpretation = TRUE)
#> When this person has a 1 probability of answering an item correctly, a same-age peer of average ability has a 0 probability of answering it correctly.
```
