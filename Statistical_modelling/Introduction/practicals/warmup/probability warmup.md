### Inheriting DNA

**Question.** Suppose <em>v<sub>1</sub></em> and <em>v<sub>2</sub></em> are two genetic variants on
the same chromosome and that there is 30% chance of recombination between them in each meiosis.
Sofiane and Ignacy are siblings. What is the probability they have inherited the same copy of
<em>v<sub>1</sub></em> and <em>v<sub>2</sub></em> (and the segment in between) from their mother?

**Note.** 'Identity by descent' (or **IBD**) is the name used for DNA segments that are identical because they are inherited from the same copy in this way.

<img src="https://render.githubusercontent.com/render/math?math=P(IBD) = P(v_1 \text{inherited IBD, and no recombination}) = P(\text{no recombination}|v_1 \text{inherited IBD}) \cdot P(\text{IBD at $v1$}) = P(\text{no recombination}) \cdot 0.5 = 0.7 * 0.5 = 0.35">

### Testing COVID

**Question.** Covid-19 lateral flow tests currently in use [are thought to be pretty
accurate](https://www.ox.ac.uk/news/2020-11-11-oxford-university-and-phe-confirm-lateral-flow-tests-
 show-high-specificity-and-are). According to that article, the false positive rate is around 0.32%:

<img src="https://render.githubusercontent.com/render/math?math=P(\text{positive}|\text{not infected}) = 0.0032">

and the specificity is around 80%:

<img src="https://render.githubusercontent.com/render/math?math=P(\text{positive}|\text{infected}) = 0.8">

Currently [about 0.4% of people in Oxfordshire are infected](https://phdashboard.oxfordshire.gov.uk):

<img src="https://render.githubusercontent.com/render/math?math=P(\text{infected}) = 0.004">

Suppose you test positive.  How likely are you to have COVID?

**Hint.** you will need Bayes' rule.
