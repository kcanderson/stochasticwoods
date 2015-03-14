(ns stochasticwoods.geneexp
  (:use [stochasticwoods.core]
        [incanter.core]
        [incanter.charts]
        [incanter.stats]))

(defn square
  [x]
  (* x x))

(defn var-imp-over-runs
  ([all_imps var]
   (map var all_imps)))

(defn var-rank-over-runs
  ([all_imps var]
   (let [imps (var-imp-over-runs all_imps var)]
     (map #(let [foo (zipmap (map first (sort-by second %1)) (range))]
             (if (zero? %2) 0
                 (foo var)))
          all_imps imps))))

(defn grab-top-n-over-runs
  [n all_imps]
  (mapcat #(map first (take n (sort-by second > %))) all_imps))

(defn grab-auc-over-runs
  [data all_forests category]
  (pmap #(classification-auc
          (oob-predictions % data)
          (map second data) category)
       all_forests))
(defn get-best-vimp
  [all_imps var]
  (loop [[imps & rimps] all_imps
         best 0]
    (cond
      (nil? imps) best
      (zero? (get imps var 0)) (recur rimps best)
      :else (recur rimps (get imps var best)))))

(defn mean-vimp
  [all_imps var]
  (mean (map #(get % var 0) all_imps)))

(defn last-usage
  [genelist included]
  (let [m (zipmap genelist (range))]
    (fn [k]
      (dec (count (remove false?
                          (map #(contains? % (m k)) included)))))))

(defn smooth-var
  [all_imps included genelist]
  (let [;;m (zipmap (range) gene_list)
        ;;k (m var)
        ;;m (zipmap gene_list (range))
        n (count included)
        luf (memoize (last-usage genelist included))]
    (fn [k]
      (let [l (luf k)
            c (map #(get % k 0) all_imps)        
            f (take (inc l) c)]
        (concat f (repeat (- n (inc l)) (mean c)))))))

(defn area-vimp
  [all_imps var]
  (area-under-curve
   (map #(list %1 (get %2 var 0)) (range) all_imps)))

(defn smooth-area-vimp
  [all_imps included]
  (let [svf (smooth-var all_imps included)]
    (fn [var]
      (area-under-curve
       (map list (range) (svf var))
       ;;(map #(list %1 (get %2 var)) (range) all_imps)
       ))))

(defn relative-ranks
  [all_imps vars_used genelist]
  (let [update (fn [imps vars]
                 (let [z (zipmap (range) genelist)
                       foo (apply merge
                                  (map #(hash-map (z %) (get imps (z %)))
                                       vars))
                       n (count foo)]
                   (reduce (fn [m [[k v] i]]
                             (assoc m k (float (/ i n))))
                           foo
                           (map list
                                (sort-by second > foo)
                                (range)))))
        ranks (map update all_imps vars_used)]
    (fn [var]
      (map #(get % var 0) ranks))))

(defn absolute-ranks
  [all_imps vars_used genelist]
  (let [update (fn [imps vars]
                 (let [z (zipmap (range) genelist)
                       foo (apply merge
                                  (map #(hash-map (z %) (get imps (z %)))
                                       vars))
                       n (count foo)
                       ;;n2 (- (count imps) n)
                       ]
                   (reduce (fn [m [[k v] i]]
                             (assoc m k i))
                           foo
                           (map list
                                (sort-by second > foo)
                                (range)))))
        ranks (map update all_imps vars_used)]
    (fn [var]
      (map #(get % var 0) ranks))))

(defn relative-rank-area-vimp
  [all_imps vars_used]
  (let [f (relative-ranks all_imps vars_used)]
    (fn [var]
      (area-under-curve
       (map #(list %1 %2) (range) (f var))))))

(defn rank-area-vimp
  [all_imps vars_used rank_fn]
  (let [f (rank_fn all_imps vars_used)]
    (fn [var]
      (area-under-curve
       (map #(list %1 %2) (range) (f var))))))

(defn special-read-string
  [x]
  (let [v (read-string x)]
    (if (symbol? v) (keyword v)
        v)))

;; Open some data
(def f "/home/kc/Code/Bioinformatics/lung_cancer/data/rnaseq_recurrence.csv")
(def f "/home/kc/Code/Bioinformatics/lung_cancer/data/lusc_recurrence_rnaseq.csv")
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/rnaseq_recurrence.csv")
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/lusc_recurrence_rnaseq.csv")
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/luad_metastasis_rnaseq.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 meta (mapv keyword (clojure.string/split (first l) #", "))
                 f (fn [row] (mapv special-read-string
                                  (clojure.string/split row #", ")))]
             [meta (doall (mapv f (rest l)))])))
(def gene_list (into [] (rest (first foo))))
(def data (map #(into [] (rest %)) (second foo)))
(def response (map first (second foo)))
(def response (map #(cond (= % :M1a) :M1
                          (= % :M1b) :M1
                          :else %)
                   response))

(def foo nil)
(def dd (map list data response))

(def f "/home/kanderson/Code/Bioinformatics/lung_risk/missing_filled_with_year.csv")
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/lusc_clinical_missing_filled.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 ks (map keyword
                         (clojure.string/split (first l)
                                               #","))
                 f (fn [row] (map special-read-string
                                 (clojure.string/split row #",")))
                 f2 (fn [row] (zipmap ks (f row)))]
             (mapv f2 (rest l)))))
(def data2 (mapv #(dissoc % :new_tumor_event_after_initial_treatment)
                 foo))
(def response2 (mapv :new_tumor_event_after_initial_treatment foo))
(def features (keys (first data2)))
(def data2 (map #(mapv (fn [k] (get % k))
                       features)
                data2))
(def dd2 (map list data2 response2))

(def f "/home/kanderson/Code/Bioinformatics/lung_risk/luad_clinical_metastasis.csv")
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/lusc_clinical_metastasis.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 ks (map keyword
                         (clojure.string/split (first l)
                                               #","))
                 f (fn [row] (map special-read-string
                                 (clojure.string/split row #",")))
                 f2 (fn [row] (zipmap ks (f row)))]
             (mapv f2 (rest l)))))
(def datam (mapv #(dissoc % :pathologic_M)
                 foo))
(def responsem (mapv #(cond
                        (= :M1a (:pathologic_M %)) :M1
                        (= :M1b (:pathologic_M %)) :M1
                        :else (:pathologic_M %))
                     foo))
(def featuresm (keys (first datam)))
(def datam (map #(mapv (fn [k] (get % k))
                       featuresm)
                datam))
(def ddm (map list datam responsem))
(def forestm (grow-forest ddm))

(def f "/home/kanderson/Downloads/glop.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 meta (mapv keyword (clojure.string/split (first l) #", "))
                 f (fn [row] (mapv special-read-string
                                  (clojure.string/split row #", ")))]
             [meta (doall (mapv f (rest l)))])))

(def gene_list (into [] (rest (first foo))))
(def data (map #(into [] (rest %)) (second foo)))
(def response (map first (second foo)))
(def foo nil)
(def dd (map list data response))

(def m0classification (grab-auc-over-runs dd allforests :M0))
(def m1classification (grab-auc-over-runs dd allforests :M1))
(def mxclassification (grab-auc-over-runs dd allforests :MX))

(use 'incanter.pdf)
(view
 (add-lines
  (add-lines (xy-plot (range 1 (inc (count allforests)))
                      m0classification
                      :series-label "M0"
                      :x-label "Runs"
                      :y-label "Classification performance (AUC)"
                      :legend true)
             (range 1 (inc (count allforests)))
             m1classification
             :series-label "M1")
  (range 1 (inc (count allforests)))
  mxclassification
  :series-label "MX"))

;; Var imp
(let [mostimps (map first (sort-by second >
                                   (filter #(> (second %) 0)
                                           (last (butlast (butlast (butlast allimps)))))))]
  (let [p (xy-plot (range 1 (inc (count allimps)))
                   (map #(/ % 1000)
                        (var-imp-over-runs allimps (first mostimps)))
                   :series-label (first mostimps)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "VIMP"                   
                   :title "Varible importance over 12 stepwise runs (top 78 genes)")]
    (view (reduce #(add-lines %1
                              (range 1 (inc (count allimps)))
                              (map (fn [i] (/ i 1000))
                                   (var-imp-over-runs allimps %2))
                              :series-label %2)
                  p (rest mostimps)))))

;; Ranks
(let [mostimps (map first
                    (sort-by second >
                             (filter #(> (second %) 0)
                                     (last (butlast (butlast allimps))))))
      ]
  (let [p (xy-plot (range (count allimps))
                   (var-rank-over-runs allimps (first mostimps))
                   :series-label (first mostimps)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "Importance rank"                   
                   :title "Varible rank over 12 stepwise runs (top 38 genes)")]
    (view (reduce #(add-lines %1
                              (range (count allimps))
                              (var-rank-over-runs allimps %2)
                              :series-label %2)
                  p (rest mostimps)))))


(def areas
  (let [genemap (zipmap (range) gene_list)]
    (map #(area-under-curve
           (map list (range)
                (var-imp-over-runs allimps (genemap %))))
     (range (count genemap)))))

(view (histogram areas))
(view (xy-plot (range (count gene_list))
               (map second
                    (sort-by second > (zipmap gene_list areas)))))

(def area_forest
  (let [a (filter #(> (second %) 7.5)
                  (zipmap (range) areas))
        features (map first a)]
    (grow-forest dd :var_sample_fn (var-sample-fn features
                                                  (java.lang.Math/sqrt
                                                   (count features))))))
;; (defn- random-split-all-vars
;;   ([data freqs split_fns gini]
;;    (random-split-all-vars data freqs split_fns
;;                           (range (count (first (first data))))))
;;   ([data freqs split_fns gini vars]
;;    (let [l (map #(cons % ((nth split_fns %) data freqs %))
;;                 vars)
;;          lf (filter #(< (second %) gini) l)
;;          n (count lf)]
;;      (if (zero? n)
;;        (list 0 1 1 1)
;;        (nth lf (rand-int n))))))

;; (defn extra-random-decision-tree
;;   ([data]
;;    (extra-random-decision-tree
;;     data (var-sample-fn data)))
;;   ([data vars_sel_fn]
;;    (extra-random-decision-tree
;;     data vars_sel_fn (gini-impurity (map second data))
;;     (build-split-fns data)))
;;   ([data vars_sel_fn score split_fns]
;;    (if (zero? score)
;;      (frequencies (map second data))
;;      (let [freqs (frequencies (map second data))
;;            [var gi gil gig
;;             lte_data gt_data
;;             split_val] (random-split-all-vars
;;                         data freqs split_fns score
;;                         (vars_sel_fn))]
;;        (if (or (and (<= score gil) (<= score gig) (<= score gi))
;;               (close-enough? gi score))
;;          freqs
;;          (list (list var split_val)
;;                (extra-random-decision-tree
;;                 lte_data vars_sel_fn gil split_fns)
;;                (extra-random-decision-tree
;;                 gt_data vars_sel_fn gig split_fns)))))))

(defn collect-stepwise-data
  ([data]
   (collect-stepwise-data
    data (range (count (first (first data)))) 0.8))
  ([data initial_variables keep_fraction]
   (loop [vars initial_variables
          forests []
          imps []
          included [(into #{} vars)]]
     (if (< (count vars) 10)
       (list forests imps included)
       (let [n (count vars)
             rf (grow-forest
                 data
                 :ntrees (max 300 (int (/ (count vars) 5)))
                 ;;:decision_tree_fn extra-random-decision-tree
                 :var_sample_fn (var-sample-fn
                                 vars (int (java.lang.Math/sqrt n))))
             ;;foo (println (count rf))
             imp (all-variables-importance rf data vars)
             ;;nnext (int (/ n 2))
             nnext (int (* keep_fraction n))
             varsnext (take nnext
                            (map first
                                 (sort-by second > imp)))]
         (recur varsnext
                (conj forests rf)
                (conj imps imp)
                (conj included (into #{} varsnext))))))))

;; (defn close-or-less-than
;;   [a b]
;;   (< a (* 1.025 b)))

;; (defn collect-stepwise-data
;;   [data]
;;   (loop [vars (range (count (first (first data))))
;;          forests []
;;          imps []
;;          included [(into #{} vars)]]
;;     (let [n (count vars)
;;           rf (grow-forest
;;               data
;;               :ntrees (max 500 (int (/ (count vars) 8)))
;;               :var_sample_fn (var-sample-fn
;;                               data vars (int (java.lang.Math/sqrt n))))
;;           imp (all-variables-importance rf data)
;;           varsnext (map first
;;                         (filter #(> (second %) 0)
;;                                 (zipmap (range) imp)))]
;;       (println (count vars) (count varsnext))
;;       (if (close-or-less-than (count vars) (count varsnext))
;;         (list forests imps included)
;;         (recur varsnext
;;                (conj forests rf)
;;                (conj imps imp)
;;                (conj included (into #{} varsnext)))))))

(do
  (def alldata (collect-stepwise-data dd))
  (def alldata2 (collect-stepwise-data dd)))

(do
  (def alldata3 (collect-stepwise-data dd))
  (def alldata4 (collect-stepwise-data dd)))

(def imps (clojure.set/union
           (last (butlast (last alldata)))
           (last (butlast (last alldata2)))
           (last (butlast (last alldata3)))
           (last (butlast (last alldata4)))))

(def imps (clojure.set/union (last (butlast (last alldata)))
                             (last (butlast (last alldata2)))))


(def allimps (map #(zipmap (map first %) (map second %)) (second alldata )))
(def allimps2 (map #(zipmap (map first %) (map second %)) (second alldata2)))

(def sds1 (map #(sd (map second %)) allimps))
(def sds2 (map #(sd (map second %)) allimps2))
(def fallimps (filter #(< (sd (map second %)) (* 7 (apply min sds1)))
                      allimps))
(def fallimps2 (filter #(< (sd (map second %)) (* 7 (apply min sds2)))
                       allimps2))

(count (clojure.set/intersection
        (into #{} (map first (filter #(> (second %) 0) (last allimps))))
        (into #{} (map first (filter #(> (second %) 0) (last allimps2))))))
(spit "lusc_rnaseq_recurrence_12runs.edn" (pr-str alldata))
(def alldata (read-string (slurp "lusc_rnaseq_recurrence_12runs.edn")))
(def alldata2 (read-string (slurp "lusc_rnaseq_recurrence_12runs2.edn")))
(def alldata2 (collect-stepwise-data dd))

(let [mostimps (map first
                    (sort-by second >
                             (filter #(> (second %) 0)
                                     (last (butlast (butlast allimps))))))]
  (let [p (xy-plot (range (count allimps))
                   (var-rank-over-runs allimps (first mostimps))
                   :series-label (first mostimps)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "Importance rank"                   
                   :title "Varible rank over 12 stepwise runs (top 40 genes)")]
    (view (reduce #(add-lines %1
                              (range (count allimps))
                              (var-rank-over-runs allimps %2)
                              :series-label %2)
                  p (rest mostimps)))))

(let [mostimps (map first (sort-by second >
                                   (filter #(> (second %) 0)
                                           (last (butlast (butlast (butlast fallimps)))))))]
  (let [p (xy-plot (range 1 (inc (count fallimps)))
                   (map #(/ % 1000)
                        (var-imp-over-runs fallimps (first mostimps)))
                   :series-label (first mostimps)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "VIMP"                   
                   :title "Varible importance over 12 stepwise runs (top 80 genes)")]
    (view (reduce #(add-lines %1
                              (range 1 (inc (count fallimps)))
                              (map (fn [i] (/ i 1000))
                                   (var-imp-over-runs fallimps %2))
                              :series-label %2)
                  p (rest mostimps)))))

(def allforests (first alldata))
(def recurrence_classification (grab-auc-over-runs dd allforests :YES))
(view (xy-plot (range 1 (inc (count allforests)))
               recurrence_classification
               :x-label "Runs"
               :y-label "Classification performance (AUC)"
               :title "LUAD recurrence stepwise performance"))




(def keylist
  (map first
       (take 20
             (sort-by second >
                      (map #(list % (area-vimp allimps %))
                           (keys (first allimps)))))))
(def keylist (keys (first allimps)))

;; Normal RF correlation
(let [k keylist]
  (view
   (scatter-plot
    (map #(get (first allimps) % 0) k)
    (map #(get (first allimps2) % 0) k)
    :x-label "Breiman VIMP (run 1)"
    :y-label "Breiman VIMP (run 2)"
    :title "LUAD VIMP correlation for two runs r2=0.006")))

(let [k keylist]
  (square
   (correlation
    (map #(/ (get (first allimps) % 0) 1000.0) k)
    (map #(/ (get (first allimps2) % 0) 1000.0) k))))


;; "Best" VIMP correlation
(let [k keylist]
  (view
   (scatter-plot
    (map #(get-best-vimp allimps %) k)
    (map #(get-best-vimp allimps2 %) k)
    :x-label "VIMP"
    :y-label "VIMP"
    :title "Last VIMP score correlation for two runs r2=0.002")))

;; Mean VIMP score
(let [k keylist]
  (println
   (square
    (incanter.stats/correlation
     (map (partial mean-vimp fallimps) k)
     (map (partial mean-vimp fallimps2) k)))))

(let [k keylist]
  (view
   (scatter-plot
    (map (partial mean-vimp fallimps) k)
    (map (partial mean-vimp fallimps2) k)
    :x-label "VIMP"
    :y-label "VIMP"
    :title "Mean VIMP score correlation for two runs r2=0.034")))

;; Area VIMP score
(let [k (keys (first allimps))]
  (println
   (square
    (incanter.stats/correlation
     (map (partial area-vimp allimps) k)
     (map (partial area-vimp allimps2) k)))))
(defn div1000
  [x]
  (/ x 1000.0))
(let [k (keys (first allimps))]
  (view
   (scatter-plot
    (map (partial area-vimp (nthrest allimps 10)) k)
    (map (partial area-vimp (nthrest allimps2 10)) k))))

(defn sig-gene?
  [[x y]]
  (and (> x 0.06) (> y 0.06)))

(defn stringify-gene-keyword
  [k]
  (first (clojure.string/split
          (subs (str k) 1) #"\|")))

(let [k keylist
      r1 (map #(/ % 1000) (map (partial area-vimp allimps) k))
      r2 (map #(/ % 1000) (map (partial area-vimp allimps2) k))]
  (def glist (map ;;(comp stringify-gene-keyword first)
              first
              (filter #(sig-gene? (rest %))
                      (map list k r1 r2)))))
(let [k keylist
      f1 (partial area-vimp allimps)
      f2 (partial area-vimp allimps2)
      r1 (map f1 k) ;;(map #(/ % 1000) (map f1 k))
      r2 (map f2 k) ;;(map #(/ % 1000) (map f2 k))
      p (scatter-plot
         r1 r2
         :x-label "Area VIMP measure (run 1)"
         :y-label "Area VIMP measure (run 2)"
         :title "LUAD area VIMP correlation for two runs r2=0.5")]
  (view (reduce (fn [c [x y i]]
                  (if (sig-gene? [x y])
                    (add-pointer c x y :text (stringify-gene-keyword (gene_list i)))
                    c))
                p (map list r1 r2 k))))

(let [k keylist]
  (view
   (add-points
    (add-points
     (scatter-plot
      (map #(/ % 1000.0) (map (partial area-vimp fallimps) k))     
      (map #(/ % 1000.0) (map (partial area-vimp fallimps2) k))
      :x-label "Area VIMP measure (run 1)"
      :y-label "Area VIMP measure (run 2)"
      :title "Area VIMP score correlation for two runs (after removal of high-variance runs)")
     (map #(/ % 1000.0) (map (partial area-vimp fallimps) cunion))     
     (map #(/ % 1000.0) (map (partial area-vimp fallimps2) cunion)))
    (map #(/ % 1000.0) (map (partial area-vimp fallimps) glist))     
    (map #(/ % 1000.0) (map (partial area-vimp fallimps2) glist)))))


(let [k keylist
      f1 (smooth-area-vimp allimps
                           (last alldata))
      f2 (smooth-area-vimp allimps2
                           (last alldata2))
      r1 (map #(/ % 1000) (map f1 k))
      r2 (map #(/ % 1000) (map f2 k))
      p (scatter-plot
         r1 r2
         :x-label "Area VIMP measure (run 1)"
         :y-label "Area VIMP measure (run 2)"
         :title "LUSC smoothed area VIMP score correlation for two runs r2=0.35")]
  (view (reduce (fn [c [x y i]]
                  (if (sig-gene? [x y])
                    (add-pointer c x y :text (stringify-gene-keyword i))
                    c))
                p (map list r1 r2 k))))

(let [k keylist
      f1 (smooth-area-vimp allimps
                           (last alldata))
      f2 (smooth-area-vimp allimps2
                           (last alldata2))
      r1 (map #(/ % 1000) (map f1 k))
      r2 (map #(/ % 1000) (map f2 k))]
  (def glist (map (comp stringify-gene-keyword first)
                  (filter #(sig-gene? (rest %))
                          (map list k r1 r2)))))



(def pcomp (principal-components
            (matrix (grab-indices-from-data dd foovars))))
(def components (:rotation pcomp))
(def pc1 (sel components :cols 0))
(def pc2 (sel components :cols 1))
(def m (matrix (grab-indices-from-data dd foovars)))
(def r (grab-indices-from-data (filter #(= :YES (last %)) dd)
                               foovars))
(def n (grab-indices-from-data (filter #(= :NO (last %)) dd)
                               foovars))

(def x1r (mmult r pc1))
(def x1n (mmult n pc1))
(def x2r (mmult r pc2))
(def x2n (mmult n pc2))

(view (add-points
       (scatter-plot x1r x2r
                     ;;:group-by species
                     :x-label "PC1" 
                     :y-label "PC2" 
                     :title ""
                     :series-label "recurrence"
                     :legend true)
       x1n x2n :series-label "no recurrence"))

;; Relative rank areas
(let [k keylist
      f3 (rank-area-vimp allimps (third alldata) relative-ranks)
      f4 (rank-area-vimp allimps2 (third alldata2) relative-ranks)]
  (view
   (scatter-plot
    (map f3 k)
    (map f4 k)
    :x-label "VIMP"
    :y-label "VIMP"
    :title "Area VIMP score correlation for two runs r2=0.042")))

;; Absolute ranks
(let [k keylist
      f3 (rank-area-vimp allimps (third alldata) absolute-ranks)
      f4 (rank-area-vimp allimps2 (third alldata2) absolute-ranks)]
  (view
   (scatter-plot
    (map f3 k)
    (map f4 k)
    :x-label "VIMP"
    :y-label "VIMP"
    :title "Area VIMP score correlation for two runs r2=0.37")))
(let [k keylist
      f3 (rank-area-vimp allimps (third alldata) absolute-ranks)      
      f4 (rank-area-vimp allimps2 (third alldata2) absolute-ranks)]
  (square (correlation
           (map f3 k)    
           (map f4 k))))

(let [k keylist
      f3 (rank-area-vimp allimps (third alldata) absolute-ranks)      
      f4 (rank-area-vimp allimps2 (third alldata2) absolute-ranks)]
  (view (histogram
         (map f4 k)
         ;;(map #(* 0.5 (+ %1 %2)) (map f3 k) (map f4 k))
           ;;(map f4 k)
           )))

(let [k (map (zipmap (range) gene_list) (last (butlast (third alldata))))
      f3 (absolute-ranks allimps3 (third alldata))
      ;;f4 (absolute-ranks allimps4 (third alldata2))
      ]
  (let [p (xy-plot (range (count allimps3))
                   (f3 (first k))
                   :series-label (first k)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "Importance rank"                   
                   :title "Varible rank over 12 stepwise runs (top 38 genes)")]
    (view (reduce #(add-lines %1
                              (range (count allimps))
                              (f3 %2)
                              :series-label %2)
                  p (rest k)))))


(let [m (zipmap (range) gene_list)
      mostimps (map m (last (third alldata)))
      f (relative-ranks allimps m)]
  (let [p (xy-plot (range 1 (inc (count allimps)))
                   (f (first mostimps))
                   :series-label (first mostimps)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "VIMP"                   
                   :title "Varible importance over 12 stepwise runs (top 80 genes)")]
    (view (reduce #(add-lines %1
                              (range 1 (inc (count allimps)))
                              (f %2)
                              :series-label %2)
                  p (rest mostimps)))))


(defn correlation-between-genes
  [data g1]
  (let [m (zipmap gene_list (range))
        d1 (map first (grab-indices-from-data data [(m g1)]))]
    (fn [g2]
      (let [d2 (map first (grab-indices-from-data data [(m g2)]))]
        (correlation d1 d2)))))
(defn correlation-between-vars
  [data v1]
  (let [d1 (map first (grab-indices-from-data data [v1]))]
    (fn [v2]
      (let [d2 (map first (grab-indices-from-data data [v2]))]
        (correlation d1 d2)))))

(defn find-correlated-vars
  [data all_vars correlation_vars threshold]
  (apply merge
         ;;(remove (comp empty? second first))
         (map (fn [g1]
                (let [f (correlation-between-vars data g1)]
                  {g1 (filter #(and (not= g1 (first %))
                                  (> (square (second %)) threshold))
                              (pmap (fn [g2] [g2 (f g2)])
                                    all_vars))}))
              correlation_vars)))
(defn union-of-correlated-vars
  [cvars]
  (into #{}
        (mapcat (fn [[k l]]
                  (map first l))
                cvars)))
(defn correlated-to-vars?
  [data var_list curr_var threshold]
  (not (empty?
      (filter #(> (square %) threshold)
              (map (correlation-between-vars data curr_var)
                   var_list)))))
;; 
(defn glopp
  [data ranked_vars corr_threshold]
  (loop [included (take 1 ranked_vars)
         [curr & remaining] (rest ranked_vars)]
    (cond
      (nil? curr) included
      (correlated-to-vars?
       data included curr corr_threshold) (recur included remaining)
       :else (recur (conj included curr) remaining)
       )))

;; 1. Find correlated vars.
;; 2. If not enough, end.
;; 3. Stepwise RF with correlated vars.
;; 4. Return to 1.
(defn glop
  [data imp_vars all_vars]
  (loop [imps [imp_vars]
         remaining_vars (clojure.set/difference all_vars imp_vars)
         correlations []]
    (let [i (last imps)
          corrs (find-correlated-vars data remaining_vars i 0.5)
          correlated (union-of-correlated-vars corrs)]
      (if (< (count correlated) 15)
        [imps correlations]
        (let [[f v vv] (collect-stepwise-data data correlated 0.75)
              perf (first (grab-auc-over-runs data [(last f)] :YES))]
          (println perf (count correlated) (count remaining_vars) (last vv))
          (if (< perf 0.6)
            [imps correlations]
            (recur (conj imps (last vv))
                   (clojure.set/difference remaining_vars (last vv))
                   (conj correlations corrs))))))))

(defn build-graph
  [m]
  (mapcat (fn [[g1 l]]
            (map (fn [[g2 v]] [g1 g2])
                 l))
          m))

(def foo (glop dd imps (into #{} (range (count gene_list)))))
(doseq [a (map (zipmap (range) gene_list) (apply concat (first foo)))]
  (println (stringify-gene-keyword a)))
(count (apply concat (first foo)))
(def corrlist (into #{} (apply concat (first foo))))
(def edges
  (filter (fn [[g1 g2]] (and (contains? corrlist g1)
                          (contains? corrlist g2)))
          (mapcat build-graph (second foo))))

(spit "/home/kc/lusc.sif"
 (with-out-str
   (let [f (comp stringify-gene-keyword (zipmap (range) gene_list))]
     (doseq [[g1 g2] edges]
       (println (f g1) " coexpressed_with " (f g2))))))

(def imps (map first imps_salient))
(def corrvars
  (find-correlated-vars dd (range (count gene_list)) imps 0.5))

(def corrs
  (apply merge
         (map (fn [g1]
                (let [f (correlation-between-genes dd g1)]
                  {g1 (filter #(and (not= g1 (first %))
                                  (> (square (second %)) 0.5))
                              (map (fn [g2] [g2 (f g2)])
                                   gene_list))}))
              (map (zipmap (range) gene_list) imps))))

(def corrimps
  (into #{} (mapcat (fn [[k l]] (map first l)) corrs)))

(count vars)
(def vars (clojure.set/union (union-of-correlated-vars corrvars) imps))
(def method1data2 (collect-stepwise-data dd imps 0.9))
(def method1data (collect-stepwise-data dd vars 0.9))
(def method1aucs (grab-auc-over-runs dd (first method1data) :YES))
(println (grab-auc-over-runs dd (first method1data2) :YES))

(view (xy-plot (range (count aucs)) aucs
               :x-label "Iterations"
               :y-label "AUC"
               :title "Stepwise RF with correlated vars"))


(classification-auc (oob-predictions (last (first alldata)) dd)
                    (map second dd)
                    :YES)
(classification-auc (oob-predictions rf2 dd)
                    (map second dd)
                    :YES)

(def allvars
  (clojure.set/union imps
                     (map (zipmap gene_list (range)) corrimps)))
(def randrf (grow-forest dd
                         :var_sample_fn (var-sample-fn
                                         (sample-without-replacement (range (count (first (first dd)))) 69)
                                         (int (java.lang.Math/sqrt 69)))))

(def subforest
  (grow-forest dd
               :var_sample_fn
               (var-sample-fn
                (apply concat (first foo))
                (java.lang.Math/sqrt (count (apply concat (first foo)))))))

(def steprand
  (collect-stepwise-data
   dd (sample-without-replacement (range (count (first (first dd)))) 69)))
(grab-auc-over-runs dd (first steprand) :YES)

(def steprf2 (collect-stepwise-data dd (map (zipmap gene_list (range)) corrimps)))
(grab-auc-over-runs dd (first steprf2) :YES)

;; Random subset of 70 vars: AUC = 0.55
(def rf (collect-stepwise-data dd allvars 0.8))
(println (grab-auc-over-runs dd (first alldata) :YES))
(classification-auc
 (oob-predictions (nth (first alldata) 29) dd)
 (map second dd) :YES)

(classification-auc (oob-predictions
                     (nth (first alldata) 29) dd)
                    (map second dd) :YES)

(grab-auc-over-runs dd (first alldata) :YES)
(classification-auc (oob-predictions
                     rfplus
                     dd)
                    (map second dd)
                    :YES)
(classification-auc (oob-predictions
                     rf
                     dd)
                    (map second dd)
                    :YES)

(float ((minimum-depth (last (first alldata))) 4427))
(mean (apply concat (minimum-depth-interaction (last (first alldata))
                                               3226 4427)))

(tree-minimum-depth-interaction (:tree (first (last (first alldata))))
                                3226 4427)
(tree-minimum-depth-interaction (:tree (first (last (first alldata))))
                                3226 3427)
(tree-minimum-depth-interaction '((:a 1) ((:b 2) {} {}) ((:b 2) {} {}))
                                :a :b)
(println (:tree (first (last (first alldata)))))
(doseq [a (clojure.set/union
           (into #{} (map (zipmap (range) gene_list) imps))
           corrimps)]
  (println (stringify-gene-keyword a)))

(doseq [a (clojure.set/union
           (into #{} (map (zipmap (range) gene_list) imps))
           ;;corrimps
           )]
  (println (stringify-gene-keyword a)))


(def bar (minimum-depth (first (first alldata))))
(minimum-depth (first (first alldata)))
(tree-minimum-depth
 (:tree (first (first (first alldata)))))

(def imps (all-variables-importance f dd))

(spit "tempstuff.edn" (pr-str [alldata alldata2]))
(def foo (read-string (slurp "tempstuff.edn")))
(def alldata (first foo))
(def alldata2 (second foo))

(def foo
  (let [ma (zipmap (range) glist)
        mb (zipmap (range) glist)]
    (fn [a b]
      ((correlation-between-genes dd (ma (int a))) (mb (int b))))))

(view (incanter.charts/heat-map foo 0 (count glist) 0 (count glist)))
(count glist)
(def rf (grow-forest dd :var_sample_fn (var-sample-fn imps 4)))
(def rf2 (grow-forest dd
                      :var_sample_fn
                      (var-sample-fn
                       (map (zipmap gene_list (range)) corrimps)
                       (int (java.lang.Math/sqrt (count corrimps))))))
(def allvars (clojure.set/union (union-of-correlated-vars corrvars)
                                imps))

(def rfplus (grow-forest dd
                         :var_sample_fn
                         (var-sample-fn
                          allvars
                          (int (java.lang.Math/sqrt (count allvars))))))


(def allimps (map #(zipmap (map first %) (map second %)) (second alldata)))
(def allimps2 (map #(zipmap (map first %) (map second %)) (second alldata2)))

(def sds1 (map #(sd (map second %)) allimps))
(def sds2 (map #(sd (map second %)) allimps2))
(def fallimps (filter #(< (sd (map second %)) (* 10 (apply min sds1)))
                      allimps))
(def fallimps2 (filter #(< (sd (map second %)) (* 10 (apply min sds1)))
                      allimps2))


(def imps_areas (map (partial area-vimp allimps) (range (count gene_list))))
(def imps_areas2 (map (partial area-vimp allimps2) (range (count gene_list))))
(def imps_mean (map #(* 0.5 (+ %1 %2)) imps_areas imps_areas2))
(def imps_salient (sort-by second >
                           (filter (fn[[i v]] (> v 0.01))
                                   (zipmap (range) imps_areas))))

(view (histogram imps_areas
                 :nbins 50))

(def cool_vars (glopp dd (map first imps_salient) 0.05))

(def cool_rf
  (grow-forest dd :ntrees 1000
               :var_sample_fn (var-sample-fn (map first imps_salient)
                                             (java.lang.Math/sqrt (count imps_salient)))))

(defn gloppy
  [thresh]
  (let [x (glopp dd (map first imps_salient) thresh)]
    [x (classification-auc
        (oob-predictions (grow-forest dd
                                      :ntrees 1000
                                      :var_sample_fn
                                      (var-sample-fn
                                       x (java.lang.Math/sqrt (count x))))
                         dd)
        (map second dd) :YES)]))

(defn glopper
  [n]
  (let [x (map first (take n imps_salient))]
    (classification-auc
     (oob-predictions (grow-forest dd
                                   :ntrees 1000
                                   :var_sample_fn
                                   (var-sample-fn
                                    x (java.lang.Math/sqrt (count x))))
                      dd)
     (map second dd) :YES)))

(println (grab-auc-over-runs dd (first alldata3) :YES))
(def x [0.01 0.02 0.03 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1 0.15 0.2 0.3 0.4 0.5])
(def y (pmap #(let [a (gloppy %)]
                (println (second a))
                a)
             x))

(view (xy-plot x (map second y)))
;; cool: (y3 y5 y7)
;; old: (y4 y6 y8)
(def y2 (pmap glopper (map (comp count first) y)))

(def y6 (pmap glopper (map (comp count first) y5)))
(def y8 (pmap glopper (map (comp count first) y7)))
(view (add-lines (xy-plot x (map second y)
                          :series-label "Top uncorrelated genes"
                          :legend true
                          :title "LUAD step-up uncorrelated method"
                          :x-label "r^2 threshold"
                          :y-label "Classification performance (AUC)")
                 x y2
                 :series-label "Top all genes"))


(classification-auc 
 (oob-predictions (grow-forest dd :ntrees 1000
                               :var_sample_fn
                               (var-sample-fn (map first imps_salient)
                                              (java.lang.Math/sqrt (count imps_salient))))
                  dd)
 (map second dd) :YES)

(classification-auc ;;(oob-predictions cool_rf dd)
 (oob-predictions (grow-forest dd :ntrees 1000
                               :var_sample_fn
                               (var-sample-fn (nth (last alldata) 29)
                                              (java.lang.Math/sqrt (count (nth (last alldata) 29)))))
                  dd)
 (map second dd) :YES)

(def cool_save cool_vars)


(def corrvars (find-correlated-vars dd (map first imps_salient) (map first imps_salient) 0.25))
(first imps_salient)
(count (union-of-correlated-vars corrvars))

(def outblah
  (with-out-str
    (let [blah (build-graph corrvars)
          m (zipmap (range) gene_list)
          f (comp stringify-gene-keyword m)
          l (clojure.set/union (into #{} (map first imps_salient))
                               (union-of-correlated-vars corrvars))
          glop (reduce clojure.set/union #{} (map #(hash-set (first %) (second %)) blah))
          bar (into #{} (map (partial apply hash-set) blah))]
      (doseq [a glop]
        (if (not (contains? l a))
          (println (f a))))
      (doseq [a bar]
        (println (f (first a)) " correlated_with " (f (second a)))))))












