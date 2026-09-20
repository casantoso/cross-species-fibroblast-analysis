# Output layout

Generated results are ignored by Git and use a consistent hierarchy.

```text
outputs/
├── preprocessing/
├── samap/
│   ├── inputs/
│   │   ├── celltype/
│   │   └── fibroblast/
│   ├── celltype/
│   │   ├── objects/
│   │   ├── tables/
│   │   ├── figures/
│   │   └── metadata/
│   └── fibroblast/
│       ├── objects/
│       ├── tables/
│       ├── figures/
│       └── metadata/
├── pyscenic/
│   ├── inputs/
│   ├── human/
│   │   ├── objects/
│   │   ├── tables/
│   │   ├── figures/
│   │   ├── logs/
│   │   └── metadata/
│   └── mouse/
│       ├── objects/
│       ├── tables/
│       ├── figures/
│       ├── logs/
│       └── metadata/
├── orthology/
│   ├── tables/
│   ├── results/
│   └── metadata/
├── enrichment/
│   ├── tables/
│   ├── figures/
│   └── metadata/
└── metadata/
```

The folders are created by the workflow helpers as needed.
