***Important information regarding the data and analysis files of Gyurkovics, Stafford, & Levita (2019)***

The files "flanker_merge.csv" and "simon_merge.csv" contain the data of all 118 participants for the flanker and the Simon tasks, respectively.
The columns code the following:
subid = the subject's ID number
subage = the participant's age entered at the start of the testing session
	(due to potential experimenter mistakes when entering these values, it is advisable to use the age information from the ppt_info file instead, see later)
gender = the participant's gender entered at the start of the testing session
	(due to potential experimenter mistakes when entering these values, it is advisable to use the gender information from the ppt_info file instead, see later)
group = the participant's age group entered at the start of the testing session
	(due to potential experimenter mistakes when entering these values, it is advisable to use the age information from the ppt_info file instead, see later) 
code = this variable was not used in analyses - it was generated in the MATLAB code to control for feature repetitions,
	it codes whether a given trial is odd-congruent, odd-incongruent, even-congruent, or even-incongruent
cong = codes the congruency of the trial (0 = congruent, 1 = incongruent)
targ = the identifier of the specific target image that was presented on that trial
resp = the key code of the response given by the participant
unknown = a constant value that was saved during data collection but was not used for anything
trial = the serial position of the given trial within a block
block = codes which block the trial comes from
Acc = the accuracy of the participant's response (0 = incorrect, 1 = correct)
RT = the response latency on that given trial
pause_time = the time the participant waited before moving on to the following trial
(For some participants there are two additional columns that can be disregarded, as they code the same information as the RT column.)

The file "sart_merge.csv" contains the SART data of 117 participants. One participant (subid 41) who completed the conflict tasks has no SART data because the testing computer malfunctioned and turned off during testing.
The data of one further participant (subid 107) had to be removed because the fire alarm went off during testing.
For the meaning of each column heading please see the explanations above. The following two variables have slightly different meanings:
code = codes the trial type of that trial (1 = Go, 2 = No Go, 3 = probe)
targ = the actual digit the participant saw on screen (on probe trials it takes the value "111")

The file "ppt_info.csv" contains demographic information, IQ scores, and the questionnaire responses of participants.
The columns code the following:
subid = the subject's ID number
mwDelib = score on the Deliberate Mind-Wandering subscale of the MWQ
mwSpont = score on the Spontaneous Mind-Wandering subscale of the MWQ
state_score = state anxiety score on the STAI
trait_score = trait anxiety score on the STAI
panas_pos1/2 = score on the positive subscale of the PANAS at time point 1/2
panas_neg1/2 = score on the negative subscale of the PANAS at time point 1/2
pub_score = self-reported puberty score for participants under 18
age_c = the participant's age, calculated by subtracting their birth date from the date of testing
IQ = the participant's IQ score
filter = codes whether the participant should be kept for analyses (0 = no, 1 = yes) 
	(for a list of reasons why certain participants were removed, see below)
sex = the participant's gender (0 = male, 1 = female)
grs = the age group of the participant (0 = Adult (reference), 1 = Early Adolescent, 2 = Mid-Adolescent, 3 = Late Adolescent)

Participants who were removed:
subid 2 - non-native speaker (all subsequent participants were native English speakers)
subid 15 - reported depression
subid 22 - reported depression and anxiety
subid 47 - reported learning difficulties
subid 48 - reported dyslexia
subid 77 - reported anxiety
subid 83 - reported living with autism spectrum disorder
subid 84 - reported dyslexia
subid 90 - reported anxiety
subid 93 - reported living with autism spectrum disorder
