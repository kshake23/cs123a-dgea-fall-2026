# Veedhi's presentation script

**Slides 5–9. Target: about 6 minutes.** Start after Kareem's handoff. Read the paragraphs; the bracketed lines are cues. Rehearse once with a timer and aim for roughly 5:45–6:15 so you are safely inside the 5–7 minute window. Say “A-G-E-R” and “S-P-P-one.”

## Slide 5 · First cohort | 0:00–1:15

Thanks, Kareem. The question we wanted to answer was whether the same changes would show up in another group of patients. I'll start with the first dataset so we have something to compare with.

[Point to the horizontal axis, then the dashed horizontal line.]

Each dot is one gene. If a dot is near zero on the horizontal axis, there was very little average difference between tumor and normal lung tissue. Dots to the right show higher expression in the tumor, and dots to the left show lower expression. Higher on the plot means a smaller adjusted p-value. The horizontal dashed line shows our cutoff after adjusting for all the genes we tested.

In this first study we had 33 matched pairs, and 6,223 of the 12,502 genes passed that cutoff. Of those, 436 are orange because they also met the size cutoff. One unit in either direction on a log-two scale means a twofold difference. So some of the gray dots above the line are still statistically significant. They just had a smaller change.

## Slide 6 · Second cohort | 1:15–2:25

GSE 7670 was the dataset I focused on. I checked the sample labels to see which ones actually belonged together. For example, a sample marked one N and a sample marked one T are the normal and tumor samples from the same patient. That gave us 27 complete pairs. The cell lines and other samples without a matching patient sample stayed out of the test.

[Point to 4,366 and 1,123.]

We tested the same 12,502 genes here. This time, 4,366 passed the significance cutoff, and 1,123 also met the size cutoff. So the second study had fewer significant genes overall but more orange dots. Those are different measures: one counts genes with strong enough statistical evidence, and the other also asks whether the difference was at least twofold. The counts alone don't tell us whether the same genes changed.

## Slide 7 · Individual matched patients | 2:25–3:35

[Trace a thin line from Normal to Tumor on each chart.]

Here we can see what those changes looked like in actual patients. On each chart, a thin line connects the normal and tumor measurements from one person. The thick line connects the group averages. The two charts have different vertical scales, so I'm focusing on the direction within each chart.

On the left, A-G-E-R was lower in the tumor in all 27 pairs. On the right, S-P-P-one was higher in the tumor in 26 of the 27 pairs. Some lines are steeper than others, so the size of the change varied across patients. For these two examples, though, the direction holds across almost every pair. We picked two clear examples, so the next plot checks how much of that agreement holds across all shared genes.

## Slide 8 · Agreement between studies | 3:35–4:55

[Point to the axes, then to the top-right and bottom-left quadrants.]

Each dot is back to being one gene. Going across, you see its estimated fold change in the first study. Going up, you see it in the second. Top right means higher in tumor in both studies. Bottom left means lower in tumor in both studies. The diagonal is where the two estimated changes would be exactly equal.

Across all 12,502 shared genes, the Spearman correlation was zero point seven eight three. Spearman compares ranks, so genes near the higher end in one study tended to be near the higher end in the other. The points still spread around the diagonal. The two studies can agree about the direction without finding exactly the same size of change.

We also checked the subset that passed the significance cutoff in both studies. Out of 3,659 genes in that subset, 3,655 had the same direction of change.

## Slide 9 · Conclusion and limits | 4:55–6:00

What I'd take away is that many tumor-associated expression changes held up in two independent groups of patients. S-P-P-one was higher and A-G-E-R was lower in both studies. That helps answer our original question about which changes we can see again in a separate cohort.

There are limits. Both studies used the same microarray platform, and we used expression values the original studies had already processed. We didn't adjust for things like stage or smoking history. We can say that expression differed between tumor and normal tissue. We can't say that these genes caused the cancer. A next step would be to see whether the same pattern holds in a group measured on a different platform.

The GitHub repo has the code, selected sample lists and figures so someone can repeat the analysis. Thanks. We're happy to take questions.

## Q&A · do not read as part of the talk

- **Why pair samples?** Expression can start at a different level in each person. For each gene, the paired test works with each person's tumor-minus-normal difference and asks if the mean of those differences differs from zero.
- **How do you know the code pairs by patient?** The manifest records each sample accession, pair ID and tissue group. A unit test shuffles the sample rows, then checks the result against known paired differences of 1, 2, 3 and 4. It gets a t statistic of about 3.873 and a two-sided p-value of about 0.0305.
- **Why were some profiles excluded?** We required exactly one tumor and one normal sample per ID. GSE10072 had 33 complete pairs and 41 excluded profiles. GSE7670 had 27 complete pairs and 12 excluded profiles, including cell lines, mixtures and commercial tissue.
- **Why use FDR?** We tested thousands of genes separately. Benjamini–Hochberg adjusts the p-values to control the expected fraction of false discoveries among the genes called significant, under the method's assumptions. The twofold size cutoff only controls the orange highlights.
- **Why average probes?** Multiple probes can map to the same gene. After removing blank or ambiguous mappings, we average retained probes per gene within each sample, then test one value per gene per sample.
- **How were the two datasets scaled?** We used the processed RMA log2 values for GSE10072 and log2-transformed the positive MAS5-scaled values for GSE7670. We analyzed the datasets separately and compared fold changes afterward.
- **Can we call the genes biomarkers?** No diagnostic or predictive performance was tested. The agreement is a check that the association reproduces across two cohorts on this platform.
