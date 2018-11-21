### 1. Assignment - Position Weight Matrices

- how many candidates are there for the different start codon variants?
- find a threshold for the logarithmic score such that 50% of the valid candidates are detected
    - parameters: length: 30, background_dist: uniform, pseudo_count: 1
    - how many false positives?
- divide data into training (first 400 sequences) and test set (rest)
- use training set to determine threshold needed to classify 50% of the valid ones right


Wie viele mögliche Kandidaten für die verschie-
denen Startcodon-Varianten gibt es insgesamt?

Versuchen Sie nun einen geeigneten Detektionsschwellwert für den logarithmi-
schen Score zu finden. Der Score wird dabei generell nur für gültige Startcodon-
Kandidaten (s.o.), vor denen mindestens L Basen beobachtbar sind, berechnet.
Wählen Sie die Schwelle t zunächst so, dass 50% der richtigen Kandidaten er-
kannt werden. Wie hoch ist der Fehler, d.h. wie viele der falschen Kandidaten
detektieren Sie damit zwangsläufig auch?


Unterteilen Sie nun die Daten in eine Trainings- und eine Test-Menge. Verwenden
Sie die ersten 400 Sequenzen zum Schätzen der Parameter und zum Ermitteln
der Detektionsschwelle gemäß dem 50% Kriterium.

Verwenden Sie den Rest der
Sequenzen, um die Erkennungsrate zu ermitteln. Wie hoch ist hier der Fehler
(falsche Kandidaten über der Schwelle)?

Mit der ROC-Kurve kann unabhängig
von einem spezifischen Schwellwert die Detektionsgenauigkeit der Methode be-
stimmt werden (s. Vorlesung). Geben Sie die Kurve aus (Funktionsplot) und
bestimmen Sie den AUC-Wert (aera under curve).