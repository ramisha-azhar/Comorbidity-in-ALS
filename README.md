# Comorbidity

## Directory structure:
```
🧠 ramisha-azhar-comorbidity-in-als/
├── 📄 README.md
├── 📁 Comorbidity.Rproj
└── 📂 code/
    ├── ⚙️ installpackages.R
    ├── 🧬 Interactome_adj.RData
    ├── 📂 files/
    │   ├── 🧾 disease_gene.xlsx
    │   ├── 📄 Phenopedia.txt
    │   └── 📄 Supplementary_data1.txt
    ├── 📊 Results/
    │   ├── 📈 matrix_separation.txt
    │   ├── 📉 matrix_separation_pval.txt
    │   ├── 📊 results_RWR.txt
    │   └── 🧠 ALS_closeness
    └── 📂 script/
        ├── 🔄 computeRWR.R
        ├── 🔗 createAdjMatrix.R
        └── 📐 module_distance_search_matrix.R
```

## **Comorbidity**

- **Definition**: 
Comorbidity is the presence of **one or more additional diseases or conditions** occurring **together with a primary disease**.

## comorbidity in ALS:  

In amyotrophic lateral sclerosis (ALS), comorbidity is an important clinical consideration, as patients frequently experience a range of neurological, psychological, respiratory, and systemic conditions that coexist with motor neuron degeneration. These comorbid conditions may arise from shared pathophysiological mechanisms, disease progression, reduced mobility, or treatment-related effects. Understanding comorbidity in ALS is essential, as it significantly influences disease management, prognosis, quality of life, and caregiver burden. Early recognition and appropriate management of comorbidities can improve symptom control, support multidisciplinary care, and optimize overall patient outcomes.

- These additional conditions may:
    - Occur **simultaneously**
    - Share **biological mechanisms**
    - Influence disease progression and treatment outcome
    - 
- Comorbidity can be:
    - **Physiological** (e.g., neurological + muscular disease)
    - **Psychological** (e.g., depression with chronic illness)

Understanding comorbidity helps in:

- Disease mechanism discovery
- Risk prediction
- Drug repurposing
  
## **Amyotrophic Lateral Sclerosis (ALS)**
- A **neurodegenerative disorder**
- Affects:
    - Upper motor neurons
    - Lower motor neurons
- Leads to:
    - Progressive muscle weakness
    - Paralysis
    - Death
- ALS disease genes:
    - Satisfy the **local and disease module hypothesis**
    - Form a **statistically significant module** in the interactome
 
## Prediction of ALS Comorbidity**

Two complementary network approaches are used:

## **1- MODULE SEPARATION ANALYSIS**

### **Module Distance Search**

- Measures how **far apart** two disease modules are in the interactome.
- Separation value = **s**

- **Interpretation**

- **s < 0**  modules overlap
- **s > 0**  modules are separated
- Statistical significance assessed using **p-value**


- **s < 0 AND p ≤ 0.05**
diseases are **significantly comorbid**

    
### **Module Distance Search Matrix**

- Computes:
    - Matrix of separation values s
    - Matrix of p-values
- Input:
    - `.xlsx` file
    - Each sheet = one disease gene list
- Output matrices are:
    - **Symmetric**
    - Compared pairwise across all diseases

---

### **ALS Module Separation Results**

ALS shows **significant overlap** with:

- Neuromuscular diseases
- Muscle weakness
- Frontotemporal dementia
- Polyneuropathies
- Vascular dementia
- Diabetic neuropathies

All have **s < 0 and p ≤ 0.05**, meaning:

- ALS module directly overlaps these disease modules
- Strong network-based evidence of comorbidity

## 2. RANDOM WALK WITH RESTART (RWR)

- Network diffusion algorithm
- Measures closeness between nodes
- Simulates a random walker moving through the interactome

### Walker behavior

At each step:
Move to a neighboring node with probability **γ**
Restart from the seed node with probability **(1 − γ)**

### Mathematical Model

<img width="303" height="72" alt="image" src="https://github.com/user-attachments/assets/7dfaca05-cefc-4972-9eea-dc071a406b97" />

Where:

- W = adjacency (transition) matrix

- E = starting vector (ALS genes = 1, others = 0)

- Rₜ = probability of being at each node at time t

- γ = restart probability (0 < γ < 1)

Iterations continue until convergence.

### Adjacency Matrix

Created using createAdjMatrix

Represents gene-gene interactions:

- 1 → interaction exists

- 0 → no interaction

Used as transition matrix **W**

Essential for RWR computation

### RWR Output

Final output = steady-state probability vector

Each gene gets a score indicating:

How close it is to ALS disease genes

### Disease Ranking Using RWR
- Step 1: Filtering

Only diseases with >30 reached genes are retained

- Step 2: ALS Closeness

For each disease:

Average RWR probabilities of its genes

Produces a disease-level scor

### Normalization

Two normalization methods:
- Z-score
  <img width="182" height="115" alt="image" src="https://github.com/user-attachments/assets/a564d12a-4ede-449c-bae4-2eb8c92102c0" />

- Modified Z-score
  <img width="343" height="132" alt="image" src="https://github.com/user-attachments/assets/8ef3db23-3155-4ecb-b00c-aa977f1f2e0c" />

Where:
- x = median
  
- MAD = median absolute deviation

Modified **z-score > 3**

Disease is considered highly associated with ALS

## text file: 
### matrix_separation.txt
<img width="1482" height="226" alt="image" src="https://github.com/user-attachments/assets/5b70f304-64fe-427c-8047-8c2b96ac314c" />

- Disease Separation Matrix

This table represents a disease separation matrix computed from the human interactome.
Each cell shows the module separation value (s) between a pair of diseases based on their disease-gene modules.

- What the Matrix Represents

Rows and columns → different diseases
(ALS, Frontotemporal Dementia, Vascular Dementia, Muscle Weakness, etc.)

Each value (s) → network-based distance between two disease modules

- The matrix is symmetric:

Distance (ALS, Muscle Weakness) = Distance (Muscle Weakness, ALS)

Diagonal values = 0

A disease compared with itself has zero separation

- How to Interpret the Values

**1. Negative values (s < 0)**

Indicate overlap or strong proximity between disease modules

Suggest shared genes, pathways, or biological mechanisms

Interpreted as potential comorbidity

- Example:

**ALS – Neuromuscular Diseases = –0.545**

Strong negative value

Indicates significant overlap

Strong evidence of comorbidity

**2. Values close to zero (s ≈ 0)**

Indicate weak or minimal overlap

Disease modules are near but mostly distinct

- Example:

**ALS – Diabetic Neuropathies ≈ –0.02**

Very small separation

Weak association

**3. Positive values (s > 0)**

Indicate separation between disease modules

Diseases are biologically distant in the interactome

Lower likelihood of comorbidity

- Example:

**Frontotemporal Dementia – Neuromuscular Diseases ≈ +0.09**

Modules are separated

Less biological overlap

- Key Observations from This Matrix

ALS shows strong overlap with:

Neuromuscular Diseases (–0.545)

Muscle Weakness (–0.321)

Frontotemporal Dementia (–0.265)

Polyneuropathies (–0.213)

These consistently negative values indicate that

ALS shares common molecular mechanisms with these diseases.

- Why This Matrix Is Important

  Quantifies disease–disease relationships

  Provides network-based evidence of comorbidity

- Helps:

 Identify related diseases

  Understand shared pathology

  Support drug repurposing studies

The disease separation matrix quantifies the network distance between disease modules. Negative separation values indicate overlapping modules and suggest potential comorbidity, while positive values indicate biologically distant diseases.

### matrix_separation_pval.txt

  
<img width="1192" height="228" alt="image" src="https://github.com/user-attachments/assets/c6fb6a7f-8ea6-4709-9ee7-90e24c7b572e" />  


“The pairwise association analysis revealed statistically significant relationships between ALS and several neuromuscular conditions, including muscle weakness, neuromuscular diseases, and polyneuropathies. A strong association was also observed between ALS and frontotemporal dementia, supporting known clinical overlap. In contrast, weaker or non-significant associations were observed between dementia subtypes and peripheral neuropathies.”

it is a pairwise statistical association matrix (most likely p-values) between ALS and several comorbid conditions. 

Rows and columns = diseases/conditions

Each cell = the statistical significance of the association between the row disease and the column disease

The values are written in scientific notation (e.g. 1.41E-17)

**- Typical interpretation**

p < 0.05  statistically significant

p ≥ 0.05  not statistically significant

## results_RWR.txt

Random Walk with Restart analysis revealed a strong network proximity between ALS and several neuromuscular and neurodegenerative disorders. As expected, motor neuron disease, hereditary spastic paraplegia, and muscular atrophy ranked among the top associated conditions, validating the network-based approach. Notably, frontotemporal dementia and frontotemporal lobar degeneration also showed high z-scores, supporting the established ALS–FTD disease spectrum. Additional associations were observed with muscle weakness, peripheral neuropathies, and vascular-related conditions, highlighting the complex and multisystem nature of ALS comorbidity.

we get the final ranked list of diseases based on their network proximity to ALS, computed using the **Random Walk with Restart (RWR) algorithm**

Each row is a disease, and the columns quantify how close that disease is to ALS in the interactome.

Diseases with the highest z-scores and modified z-scores are the strongest ALS-associated conditions, including:

Amyotrophic Lateral Sclerosis (self-validation)

Motor Neuron Disease

Hereditary Spastic Paraplegia

Muscular Atrophy

Frontotemporal Dementia

This confirms biological validity of the method.

- Disease

This column lists the diseases or clinical conditions analyzed in the study. Each disease represents a set of genes mapped onto the biological interaction network and evaluated for its relationship with ALS.

- Node_num

This column indicates the number of genes (nodes) associated with each disease in the network.
A higher node count means the disease is represented by more genes.
This provides context for the robustness and coverage of the disease within the interactome.

- Prob_mean

This represents the mean steady-state probability obtained from the Random Walk with Restart (RWR) algorithm.
It reflects how frequently the random walker visits nodes linked to a specific disease.
Higher values indicate closer network proximity to ALS.
This metric captures the overall influence of ALS seed nodes on each disease module.

- Zscore

The z-score measures how much the observed network proximity deviates from what would be expected by chance.
It is calculated by comparing the observed value to a randomized background distribution.
A z-score greater than 2 is generally considered statistically significant.
Higher z-scores indicate stronger enrichment and non-random association with ALS.

- Zscore_mod

This is a modified z-score, designed to be more robust and less sensitive to extreme values or outliers.
It provides a more stable measure for ranking diseases.
This score is particularly useful for comparative analysis across conditions.
Higher values confirm the reliability of the observed association

## ALS Closeness Plot

<img width="626" height="607" alt="image" src="https://github.com/user-attachments/assets/7b73bd84-2563-4ca7-8bb0-117e61ea09d8" />

The grey density curve represents the null distribution of ALS closeness values obtained from randomized network simulations. Colored dashed vertical lines indicate the observed closeness scores of selected diseases. Diseases located in the right tail of the distribution exhibit significantly higher network proximity to ALS than expected by chance, supporting their potential comorbidity at the molecular level.

- the distribution of ALS closeness scores obtained from the Random Walk with Restart (RWR) analysis and compares the observed disease scores to a background (null) distribution.
The grey shaded curve represents the null distribution of ALS closeness scores.

- This distribution is generated from randomized disease modules or random walks, reflecting what closeness values would be expected by chance.

- Most random values cluster around zero, indicating no meaningful network association with ALS.
  
 **X-axis: ALS closeness**

- This axis represents the network proximity score to ALS.

 - Higher values indicate stronger molecular closeness to ALS within the interactome.

- Values far to the right of the distribution suggest non-random, biologically meaningful associations.

**Y-axis: Probability density**

- This shows how frequently a given closeness value occurs in the null model.

- Peaks indicate common random outcomes, while the tail represents rare events.

The plot visually confirms that several diseases have ALS closeness scores well beyond random expectation.
This provides statistical and visual validation of the RWR-based comorbidity analysis.
Together with z-scores and p-values, it strengthens confidence that the identified diseases are true ALS-related comorbidities, not random artifacts.
