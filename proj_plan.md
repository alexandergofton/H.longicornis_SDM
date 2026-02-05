# Haemaphysalis longicornis Species Distribution Modeling Project Plan

## Phase 1: Foundation & Data Collection
### Step 1: Gather Occurrence Data

**Goal:** Get a dataset of where H. longicornis has been found in Australia
**Why:** SDMs need real locations where the species exists. Think of it like teaching a computer "this tick lives in places that look like this"
**What you'll do:** Download occurrence records from GBIF (Global Biodiversity Information Facility), filter for Australia only, clean up duplicates and dodgy coordinates
Skills you'll learn: Data downloading, basic data cleaning in R

### Step 2: Get Climate Data

**Goal:** Download current climate variables for Australia
**Why:** These are the environmental conditions that determine where the tick can survive (temperature, rainfall, etc.)
**What you'll do:** Download WorldClim bioclimatic variables (there are 19 standard ones), crop them to Australia's extent
Skills you'll learn: Working with raster data, spatial cropping

### Step 3: Data Exploration

**Goal:** Understand your data before modeling
**Why:** You need to know what you're working with - are your occurrence points clustered? Are they across diverse climates?
**What you'll do:** Make simple maps, look at climate values at occurrence points, check for spatial bias
Skills you'll learn: Basic spatial visualization, exploratory data analysis

## Phase 2: Building Your First Model
### Step 4: Variable Selection

**Goal:** Choose which climate variables to use
**Why:** Too many correlated variables confuse the model. You want variables that actually matter to tick survival
**What you'll do:** Check correlations between variables, remove highly correlated ones (like keeping Bio1 but dropping Bio5 if they're >0.8 correlated)
Skills you'll learn: Correlation analysis, ecological thinking about what matters to organisms

### Step 5: Run Your First MaxEnt Model

**Goal:** Create a basic habitat suitability map
**Why:** This is your proof of concept - can you actually make this work?
**What you'll do:** Use the 'dismo' package in R to run MaxEnt with default settings on your cleaned data
Skills you'll learn: MaxEnt basics, interpreting suitability maps

### Step 6: Evaluate Model Performance

**Goal:** Figure out if your model is any good
**Why:** A bad model is worse than no model - you need to know if it's reliable
**What you'll do:** Look at AUC scores, test-train splits, variable contribution plots
Skills you'll learn: Model evaluation, understanding what makes a good vs bad model

## Phase 3: Refinement
### Step 7: Optimize Model Settings

**Goal:** Fine-tune your model for better performance
**Why:** Default settings aren't always best. Tweaking parameters can dramatically improve accuracy
**What you'll do:** Use ENMeval package to test different feature classes and regularization settings
Skills you'll learn: Model tuning, understanding model complexity vs overfitting

### Step 8: Validate with Independent Data

**Goal:** Test if your model works on data it hasn't seen
**Why:** A model that only works on training data is useless
**What you'll do:** Hold out 20-30% of your data, see if model predictions match reality
Skills you'll learn: Cross-validation, thinking critically about model reliability

## Phase 4: Future Climate Projections
### Step 9: Get Future Climate Scenarios

**Goal:** Download climate projections for 2050 and 2070
**Why:** Your whole project goal is to see where the tick might spread under climate change
**What you'll do:** Download CMIP6 climate data for different scenarios (like SSP2-4.5 and SSP5-8.5 from WorldClim)
Skills you'll learn: Understanding climate scenarios, working with projections

## Step 10: Project Your Model to Future Climates

**Goal:** Create maps showing future suitable habitat
**Why:** This is the payoff - where will this tick be able to live in 30-50 years?
**What you'll do:** Apply your trained model to future climate layers
Skills you'll learn: Model projection, dealing with extrapolation issues

## Step 11: Visualize and Communicate Results

**Goal:** Make clear, interpretable maps and figures
**Why:** Science is only useful if you can show what you found
**What you'll do:** Create before/after maps, calculate area changes, identify high-risk regions
Skills you'll learn: Scientific visualization, storytelling with data
