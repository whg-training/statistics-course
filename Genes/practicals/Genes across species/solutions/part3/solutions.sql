-- Find extreme genes by size or number of exons or transcripts:

-- TEMP tables only exist for this connection, and will disappear as soon as we close it.
.mode column
.header on

CREATE TEMP TABLE extremes AS
SELECT analysis,
       MAX( end - start ) AS max_size,
       MIN( end - start ) AS min_size,
       MAX( average_number_of_exons ) AS max_exons,
       MAX( number_of_transcripts) AS max_transcripts
FROM gene_summary_view
GROUP BY analysis ;

-- Find genes with highest number of transcripts...
SELECT E.*, G.*
FROM extremes E
INNER JOIN gene_summary_view G
ON G.analysis == E.analysis
AND G.number_of_transcripts == E.max_transcripts ;

-- Find genes with highest length...
SELECT E.*, G.*
FROM extremes E
INNER JOIN gene_summary_view G
ON G.analysis == E.analysis
AND (G.end - G.start) == E.max_size ;

-- Find genes with highest number of exons...
SELECT E.*, G.*
FROM extremes E
INNER JOIN gene_summary_view G
ON G.analysis == E.analysis
AND G.average_number_of_exons == E.max_exons ;

-- Find proportion of single-exon genes

DROP VIEW single_exon_gene_count ;
CREATE TEMP VIEW single_exon_gene_count AS
SELECT analysis, COUNT(*) AS `count` FROM gene_summary_view WHERE biotype == 'protein_coding' AND average_number_of_exons == 1.0
GROUP BY analysis ;

DROP VIEW protein_coding_gene_count ;
CREATE TEMP VIEW protein_coding_gene_count AS
SELECT analysis, COUNT(*) AS `count` FROM gene_summary_view WHERE biotype == 'protein_coding'
GROUP BY analysis ;

CREATE TEMP TABLE analyses AS SELECT DISTINCT analysis FROM genes ;

CREATE TEMP VIEW single_exon_proportion_view AS
SELECT T.analysis, S.`count` AS single_exon_gene_count, P.`count` AS `total`, CAST(S.`count` AS FLOAT) / P.`count` AS proportion
FROM analyses T
INNER JOIN single_exon_gene_count S
ON S.analysis == T.analysis
INNER JOIN protein_coding_gene_count P
ON P.analysis == T.analysis ;

SELECT * FROM single_exon_proportion_view ;
