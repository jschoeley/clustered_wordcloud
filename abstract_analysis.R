############################################
# Co-word Analysis of MaxO Paper Abstracts #
############################################

# Biodemography is what we do. What we do is what we write. What do we write?
# This analysis shows how MaxO defines the field of Biodemography through their
# publication practice. The words in the paper abstracts and their co-occurences
# constitute a semantic network. Terms central to Biodemography will be central
# in the network while more specialised bidemographical sub-topics will
# constitute sub-networks of terms.

# Init --------------------------------------------------------------------

library(tidyverse)
library(igraph)

# Input -------------------------------------------------------------------

# maxo publication records
pub <-
  readxl::read_excel(
    "publication_list.xlsx",
    sheet = "publications",
    na = "NA"
  )
aut <-
  readxl::read_excel(
    "publication_list.xlsx",
    sheet = "authors",
    na = "NA"
  )
link <-
  readxl::read_excel(
    "publication_list.xlsx",
    sheet = "publication_author_link",
    na = "NA"
  )

linkaut <-
  left_join(link, aut, by = "author_id")

# Preparing text corpus ---------------------------------------------------

# stopwords manually determined by looking though the list of words in corpus
additional_stopwords <-
  # words I deem un-informative
  c("increase", "estimated", "lived", "compared", "related",
    "specific", "suggest", "based", "level", "significantly",
    "dependent", "environmental", "including", "observed",
    "developed", "examine", "found", "assessed", "determine",
    "result", "approach", "investigate", "understanding",
    "affect", "os", "exist", "factors", "contribute",
    "evidence", "low", "provide", "ss", "support",
    "abundance", "addition", "major", "account", "common",
    "linkage", "reduced", "revealed", "set", "term", "due",
    "reported", "total", "applications", "applied", "average",
    "discuss", "examples", "i.e", "limited", "occurred",
    "underlying", "values", "calculated", "concept",
    "considerable", "consistent", "produced", "aim", "aspect",
    "captures", "conducted", "familial", "lower", "require",
    "strong", "target", "variants", "analyzed", "approximately",
    "closely", "collected", "complete", "evaluate", "experienced",
    "expressed", "focus", "growing", "half", "illustrate",
    "interpretation", "obtained", "preferred", "previously",
    "rapidly", "reach", "status", "type", "varied", "suggested",
    "regard", "similar", "lead", "simple", "describe", "role", "e.g.",
    "e.g", "introduce", "key", "accounting", "difficult", "literature",
    "explain", "achieve", "range", "wide", "contributions",
    "identified", "established", "additional", "form", "plot",
    "directions", "combined", "recent", "examined", "reveal", "main",
    # extra stopwords found in many abstracts
    "abstract", "objective", "objectives",
    "conclusion", "conclusions", "contribution",
    "methods", "background", "results")

# preparing the text corpus
corpus <-
  pub %>%
  # turn some important 2-grams into single words
  mutate(
    publication_abstract =
      gsub(
        pattern = "life expectancy",
        replacement = "lifeexpectancy",
        publication_abstract, ignore.case = TRUE)
  ) %>%
  # split abstracts into single words (tokenization)
  # this automatically removes punctuation and whitespace
  tidytext::unnest_tokens(
    output = term, input = publication_abstract, token = "words"
  ) %>%
  # normalize words
  mutate(
    # transform to lowercase
    term =
      str_to_lower(term) %>%
      # transform special chars to standard ascii (transliteration)
      stringi::stri_trans_general("Latin-ASCII") %>%
      # remove numbers or "words" containing numbers
      {ifelse(grepl("[0-9]", .), NA, .)}
  ) %>%
  filter(!is.na(term)) %>%
  # stem the words
  mutate(
    term_stem =
      SnowballC::wordStem(words = term, language = "eng")
  ) %>%
  # replace each stemmed word with its most frequent un-stemmed form
  group_by(term_stem) %>%
  mutate(
    term_stem_proxy =
      names(sort(table(term), decreasing = TRUE))[1]
  ) %>%
  ungroup() %>%
  # remove stopwords
  anti_join(
    tidytext::stop_words,
    by = c("term_stem_proxy" = "word")
  ) %>%
  anti_join(
    tibble(word = additional_stopwords),
    by = c("term_stem_proxy" = "word")
  ) %>%
  # discard words occuring only a few times in whole corpus
  group_by(term_stem_proxy) %>%
  filter(n() > 6) %>%
  arrange(publication_id) %>%
  ungroup()
glimpse(corpus)

# merge author names into corpus
corpus <-
  corpus %>%
  group_by(publication_id) %>%
  do({
    # the maxo-authors of the given paper
    paper_maxo_aut <-
      filter(linkaut,
             publication_id == .$publication_id[1],
             grepl("*max", x = author_paper_affil)) %>%
      select(author_name_short) %>% unlist()
    tibble(
      term = c(.$term_stem_proxy, paper_maxo_aut),
      type = c(rep("topic", length(.$term_stem_proxy)),
               rep("author", length(paper_maxo_aut)))
    )
  })

# Term frequency ----------------------------------------------------------

# calculate term frequencies
term_freq <-
  corpus %>%
  group_by(term) %>%
  summarise(
    # term type
    type = type[1],
    # term frequency
    tf  = n(),
    # inverse document frequency
    idf = log(length(unique(corpus$publication_id))/length(unique(publication_id))),
    # weighted term frequency
    tf_idf = tf*idf
    ) %>%
  ungroup() %>%
  arrange(desc(tf))

# Calculate distance matrix between words ---------------------------------

# generate co-occurence matrix of terms and igraph object
term_cooc <-
  table(corpus$publication_id, corpus$term) %>%
  crossprod()

# convert co-occurence matrix into undirected graph
term_graph <-
  term_cooc %>%
  graph_from_adjacency_matrix(
    mode = "undirected", weighted = TRUE, diag = FALSE
  ) %>%
  as_data_frame() %>%
  graph_from_data_frame(directed = FALSE, vertices = term_freq) %>%
  # delete edges between authors
  # (we only want to see their relation to terms, not between each other)
  {delete.edges(
    .,
    E(.)[V(.)[type == "author"] %--% V(.)[type == "author"]]
  )}

term_distances <- distances(term_graph, algorithm = 'unweighted')

# Dimensionality reduction of distance matrix -----------------------------

# metric dimensional scaling with unweighted shortest paths matrix
# (degrees of separation) between any two words based on co-occurence relations
set.seed(4)
xy_positions <-
  layout_with_mds(
    term_graph, dim = 2,
    dist = term_distances
  ) %>% as_tibble(.name_repair = function (x) c('x', 'y'))

# Plot topically clustered wordcloud --------------------------------------

data_for_plot <-
  bind_cols(
    label = V(term_graph)$name,
    tf = V(term_graph)$tf,
    type = V(term_graph)$type,
    xy_positions
  )

# topical cluster word-cloud
# merge coocurence matrix with node attributes and make graph
ggplot(data_for_plot) +
  ggrepel::geom_text_repel(
    aes(
      x = x, y = y,
      label = label, size = tf, alpha = tf, color = type
    ),
    fontface = "bold",
    # we avoid overlap and reduce distance between nodes by
    # using a force directed layout algorithm
    force = 0.1, segment.size = NA
  ) +
  scale_size(range = c(3,12), guide = "none") +
  scale_alpha(range = c(0.2, 1), guide = "none") +
  scale_color_manual(values = c("red", "black"), guide = "none") +
  coord_fixed(ratio = 0.7) +
  theme_void()
