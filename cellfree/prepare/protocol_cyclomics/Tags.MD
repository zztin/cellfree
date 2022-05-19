# Additional SAM file information for Cyclomics

The SAM standard has been designed during the rise of Illumina sequencing, as such many of its features represent technical details coming from challenges with that technology. Due to this we make some non breaking additions.

## Additional Tags

According to the SAM spec you can add tags for local use using lowercase letters, or using any tag in the X?, Y? or Z? range. We use X for sequencing related tags. Y for consensus generation tags, and Z for other tags.

### X tags

| TAG | Type  | Description                                                                                                        | remark                          | Category       |     |
| --- | ----- | ------------------------------------------------------------------------------------------------------------------ | ------------------------------- | -------------- | --- |
| XC  | Int   | pore channel on the flowcell.                                                                                      | available in sequencing summary | Sequencer info |     |
| XM  | Int   | pore mux within channel                                                                                            | available in sequencing summary | Sequencer info |     |
| Xc  | Str   | combination of mux & channel, unique identifier.                                                                   | available in sequencing summary | Sequencer info |     |
| XT  | Float | start time (s) from start of run.                                                                                  | available in sequencing summary | Sequencer info |     |
| XD  | Float | duration (s) of the read                                                                                           | available in sequencing summary | Sequencer info |     |
| XA  | Float | adapter_duration                                                                                                   | available in sequencing summary | Sequencer info |     |
| XL  | Int   | Number of bases in the output sequence, taking into </br> account any sequence trimming (sequence_length_template) | available in sequencing summary | Sequencer info |     |
| XQ  | Float | The qscore corresponding to the mean error rate of </br>the sequence (mean_qscore_template)                        | available in sequencing summary | Sequencer info |     |
| XT  | Float | The median current of the read, in pA. (median_template)                                                           | available in sequencing summary | Sequencer info |     |
| XM  | Float | The median absolute deviation of the current of </br>the read, in pA. (mad_template)                               | available in sequencing summary | Sequencer info |     |
| XE  | Str   | reason to end sequencing                                                                                           | available in sequencing summary | Sequencer info |     |


### Y tags

| TAG | Type | Description                                           | remark | Category |     |
| --- | ---- | ----------------------------------------------------- | ------ | -------- | --- |
| YB  |      | Backbone barcode                                      |        |          |     |
| YL  |      | Length of the generated repeat                        |        |          |     |
| YC  |      | Classification of the read using the consensus caller |        |          |     |
| YT  |      | derrived original total length                        |        |          |     |
| YP  |      | location of partner                                   |        |          |     |
| YA  |      | Ancestor mapping position                             |        |          |     |
| YM  |      | (Magnitude) Number of copies in the consensus         |        |          |     |
