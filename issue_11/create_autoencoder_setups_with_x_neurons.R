model_filenames <- gcaer::get_gcae_model_filenames()

for (model_filename in model_filenames) {
  for (n_neurons in seq(1, 5)) {
    model <- gcaer::read_model_file(model_filename)
    new_model <- set_n_neurons_in_latent_layer(
      model = model,
      n_neurons = n_neurons
    )
    new_model_filename <- stringr::str_replace(
      string = basename(model_filename),
      pattern = ".json$",
      replacement = paste0("_", n_neurons, "n.json")
    )
    message(new_model_filename)
    save_model(model = new_model, model_filename = new_model_filename)
  }
}

