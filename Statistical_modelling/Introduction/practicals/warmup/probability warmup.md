## Probability warm-up questions

Can you answer (or start to answer) these using the rules on the left side of the [probability
cheatsheet](../../notes/Probability%20cheatsheet.pdf)?

### DNA inheritance

When copies of the same ancestral stretch of DNA are inherited by two individuals, we say they have
'identity by descent' (**IBD**).

**Question.** Suppose <em>v<sub>1</sub></em> and <em>v<sub>2</sub></em> are two genetic variants on
the same chromosome and that there is 30% chance of recombination between them in each meiosis.
Sofiane and Ignacy are siblings. What is the probability they have inherited the same copy of the stretch of DNA between <em>v<sub>1</sub></em> and <em>v<sub>2</sub></em> from their mother?

**Note.** When the same copy of a piece of DNA is inherited by two individuals, we say they have 'identity by descent' (**IBD**). 

**Hint.** Inheriting the whole segment IBD is the same as inheriting the DNA at the first variant IBD, and there having been no recombinations.  Use this and apply Rule III.

### COVID testing

**Question.** Covid-19 lateral flow tests currently in use [are thought to be pretty
accurate](https://www.ox.ac.uk/news/2020-11-11-oxford-university-and-phe-confirm-lateral-flow-tests-show-high-specificity-and-are). According to that article, the false positive rate is around 0.32% 
<img src="https://render.githubusercontent.com/render/math?math=P(\text{positive}|\text{not infected}) = 0.0032">

and the sensitivity is around 80%:

<img src="https://render.githubusercontent.com/render/math?math=P(\text{positive}|\text{infected}) = 0.8">

Currently [about 0.4% of people in Oxfordshire are confirmed positive every
week](https://phdashboard.oxfordshire.gov.uk). That number will only capture the people who are
tested, so let's suppose around twice that number (0.8%) are currently infected.

<img src="https://render.githubusercontent.com/render/math?math=P(\text{infected}) = 0.004">

Suppose you test positive.  How worried are you that you have COVID?

### How fair is unfair?

**Question.** I toss a coin ten times and get ten heads.  How likely is it that the coin is fair?

**Note.** The inheritance question above is of a different type than the other two. In the
inheritance question, we have a model of what's going on and are simply asked to calculate a
probability using the model. In the other two questions, we have some actual *data* - the positive
test in the first question, or the sequence of heads in the second - and we are asked to draw inference using that data. This is the usual form of of scientific questions and it is typical that Bayes' theorem is needed to answer them.
