#' @description Classes and functions for handling BIDS-formatted directories
#' @importFrom fs dir_exists file_exists dir_ls path
#' @importFrom jsonlite fromJSON
#' @importFrom RNifti readNifti
#' @importFrom purrr map_dfr map
#' @importFrom R6 R6Class

#' @title Validate path
#' @description Verifies if a path exists as a file or directory.
#' @param input_path Path to validate
#' @return TRUE upon successful validation (otherwise throws error).
#' @keywords internal
validate_path <- function(input_path) {
  if (!fs::dir_exists(input_path) && !fs::file_exists(input_path)) {
    stop(sprintf("Path not found: %s", input_path))
  }
  TRUE
}

#' @export
#' @title BIDS Directory Handler
#' @description Handles BIDS-formatted directories.
Bids <- R6::R6Class("Bids",
  private = list(
    #' @description Determine the path to the raw data folder.
    #' It checks for a 'raw' subdirectory, otherwise uses the base path.
    #' @param base_bids_path Character string of the base BIDS directory path.
    #' @keywords internal
    raw_folder = function(base_bids_path) {
      raw_path <- fs::path(base_bids_path, "raw")
      if (!fs::dir_exists(raw_path)) {
        raw_path <- base_bids_path
      }
      raw_path
    },

    #' @description Load subject information from the raw data directory.
    #' Identifies directories starting with 'sub-' and creates Subject objects.
    #' @keywords internal
    load_subjects = function() {
      subject_dirs <- fs::dir_ls(self$raw_path, type = "directory")
      subject_dirs <- subject_dirs[grepl("^sub-", basename(subject_dirs))]

      if (length(subject_dirs) == 0) {
        stop(sprintf("No subjects found in %s", self$raw_path))
      }

      self$subjects <- lapply(subject_dirs, function(dir) {
        Subject$new(dir)
      })
    },

    #' @description Load derivative data for all subjects and sessions.
    #' Looks for a 'derivatives' directory and processes each subject/session.
    #' @keywords internal
    load_derivatives = function() {
      # Use the base path to find derivatives
      derivatives_path <- fs::path(self$path, "derivatives")
      if (!fs::dir_exists(derivatives_path)) {
        return()
      }

      derivative_dirs <- fs::dir_ls(derivatives_path, type = "directory")

      for (derivative_dir in derivative_dirs) {
        for (subject in self$subjects) {
          subject_deriv_dir <- fs::path(derivative_dir, subject$get_name())
          if (!fs::dir_exists(subject_deriv_dir)) next

          for (session in subject$sessions) {
            tryCatch(
              {
                session_deriv_dir <- fs::path(
                  subject_deriv_dir,
                  session$get_name()
                )
                if (!fs::dir_exists(session_deriv_dir)) next

                session$load_scan_types(session_deriv_dir)
              },
              error = function(e) {
                warning(sprintf(
                  "Error loading derivatives for %s/%s: %s",
                  subject$get_name(), session$get_name(),
                  e$message
                ))
              }
            )
          }
        }
      }
    }
  ),
  public = list(
    #' @field path Path to the BIDS root directory.
    path = NULL,
    #' @field raw_path Path to the raw data directory.
    raw_path = NULL,
    #' @field subjects List containing subject objects found in the directory.
    subjects = list(),

    #' @description Initializes a BIDS object with the given path.
    #' @param bids_root_path Path to the BIDS directory
    initialize = function(bids_root_path) {
      validate_path(bids_root_path)
      self$path <- bids_root_path
      self$raw_path <- private$raw_folder(bids_root_path)
      private$load_subjects()
      private$load_derivatives()
    },

    #' @description Prints a tree-like structure of the BIDS directory.
    #' Shows subjects, sessions, and scans.
    #' @param include_details Logical, whether to include full file paths for
    #' scans (default: FALSE).
    print_tree = function(include_details = FALSE) {
      for (i in seq_along(self$subjects)) {
        subject <- self$subjects[[i]]
        is_last_subject <- i == length(self$subjects)
        subject_prefix <- if (is_last_subject) "└── " else "├── "
        cat(sprintf("%s%s\n", subject_prefix, subject$format()))

        session_indent <- if (is_last_subject) "    " else "│   "

        for (j in seq_along(subject$sessions)) {
          session <- subject$sessions[[j]]
          is_last_session <- j == length(subject$sessions)
          session_prefix <- if (is_last_session) "└── " else "├── "
          cat(sprintf(
            "%s%s%s\n", session_indent,
            session_prefix, session$format()
          ))

          scan_indent <- if (is_last_session) {
            paste0(session_indent, "    ")
          } else {
            paste0(session_indent, "│   ")
          }

          for (k in seq_along(session$scans)) {
            scan <- session$scans[[k]]
            is_last_scan <- k == length(session$scans)
            scan_prefix <- if (is_last_scan) "└── " else "├── "

            scan_info <- if (include_details) {
              sprintf("%s (%s)", scan$format(), scan$path)
            } else {
              scan$format()
            }

            cat(sprintf("%s%s%s\n", scan_indent, scan_prefix, scan_info))
          }
        }
      }
    },

    #' @description Set oblique center and affine transformation matrix from parcellation images matching a scan
    #' pattern for each session.
    #' @param scan_pattern Regular expression pattern to match scan names
    #' @param folder_pattern Optional regular expression pattern to match folder names
    #' @return List of Slice objects
    set_oblique = function(scan_pattern, folder_pattern = NULL) {
      for (subject in self$subjects) {
        for (session in subject$sessions) {
          matching_scans <- session$scans[grepl(scan_pattern, sapply(session$scans, function(s) s$scan_name))]
          if (!is.null(folder_pattern)) {
            matching_scans <- matching_scans[grepl(folder_pattern, sapply(matching_scans, function(s) s$path))]
          }
          if (length(matching_scans) == 0) {
            next
          }
          if (length(matching_scans) > 1) {
            warning(sprintf("Multiple parcellation images found for scan pattern: %s \nTaking first one", scan_pattern))
          }
          parcellation_img <- load_nifti(matching_scans[[1]]$path)
          oblique_affine <- compute_combined_rotation(parcellation_img)
          oblique_center <- calculate_parcellation_centroid(parcellation_img)
          session$set_oblique_affine(oblique_affine)
          session$set_oblique_center(oblique_center)
        }
      }
    },

    #' @description Get standard slices from scans matching a scan pattern
    #' @param scan_pattern Regular expression pattern to match scan names
    #' @param folder_pattern Optional regular expression pattern to match folder names
    #' @param oblique_slices Logical, whether to extract oblique slices using the oblique affine
    #' @return List of Slice objects
    get_slices = function(scan_pattern, folder_pattern = NULL, oblique_slices = FALSE) {
      result_slices <- tibble::tibble()
      for (subject in self$subjects) {
        for (session in subject$sessions) {
          matching_scans <- session$scans[grepl(scan_pattern, sapply(session$scans, function(s) s$scan_name))]
          if (!is.null(folder_pattern)) {
            matching_scans <- matching_scans[grepl(folder_pattern, sapply(matching_scans, function(s) s$path))]
          }
          if (length(matching_scans) == 0) {
            next
          }
          for (scan in matching_scans) {
            img <- load_nifti(scan$path)
            if (oblique_slices) {
              if (is.null(session$oblique_affine) || is.null(session$oblique_center)) {
                warning(sprintf("Oblique slices requested but no oblique affine or center found for session %s", session$get_name()))
                extracted_slices <- extract_nonzero_slices(img)
              }
              extracted_slices <- extract_oriented_slices(img,
                                                          affine = session$oblique_affine,
                                                          center = session$oblique_center)
            } else {
              if (is.null(session$oblique_affine) || is.null(session$oblique_center)) {
                warning(sprintf("Oblique slices requested but no oblique affine or center found for session %s", session$get_name()))
                extracted_slices <- extract_nonzero_slices(img)
              } else {
                extracted_slices <- extract_nonzero_slices(img, center = session$oblique_center)
              }
            }
            dataframe <- tibble::tibble(
              subject_id = subject$id,
              session_id = session$id,
              scan_name = scan$scan_name,
              slice_type = c("sagittal", "coronal", "axial"),
              data = list(
                extracted_slices$sagittal,
                extracted_slices$coronal,
                extracted_slices$axial
              )
            )
            result_slices <- dplyr::bind_rows(result_slices, dataframe)
          }
        }
      }

      if (nrow(result_slices) == 0) {
        stop(sprintf("No slices could be extracted for scans matching scan pattern: %s", scan_pattern))
      }

      result_slices
    }
  )
)

#' @export
#' @title BIDS Subject Handler
#' @description Represents a single subject within a BIDS dataset.
#' @keywords internal
Subject <- R6::R6Class("Subject",
  private = list(
    #' @description Load session information for the subject.
    #' Identifies directories starting with 'ses-' within the subject's path.
    #' @param subject_dir_path Character string path to the subject's directory.
    #' @return NULL. Populates the `self$sessions` list.
    #' @keywords internal
    load_sessions = function(subject_dir_path) {
      validate_path(subject_dir_path)
      session_dirs <- fs::dir_ls(subject_dir_path, type = "directory")
      session_dirs <- session_dirs[grepl("^ses-", basename(session_dirs))]

      if (length(session_dirs) == 0) {
        stop(sprintf("No sessions found in %s", subject_dir_path))
      }

      self$sessions <- lapply(session_dirs, function(dir) {
        Session$new(dir)
      })
    }
  ),
  public = list(
    #' @field path Path to the subject's directory.
    path = NULL,
    #' @field id Identifier for the subject (e.g., "01").
    id = NULL,
    #' @field sessions List containing session objects for this subject.
    sessions = list(),

    #' @description Initializes a subject object.
    #' @param subject_path Character string path to the subject's directory.
    initialize = function(subject_path) {
      validate_path(subject_path)
      self$path <- subject_path
      self$id <- sub("sub-", "", basename(subject_path))
      private$load_sessions(subject_path)
    },

    #' @description Get the subject identifier.
    #' @return Character string subject ID (e.g., "01").
    get_id = function() {
      self$id
    },

    #' @description Get the full subject name (e.g., "sub-01").
    #' @return Character string subject name with prefix.
    get_name = function() {
      sprintf("sub-%s", self$id)
    },

    #' @description Format the subject information for printing.
    #' @return Character string representation (e.g., "Subject(id=01)").
    format = function() {
      sprintf("Subject(id=%s)", self$id)
    }
  )
)

#' @export
#' @title BIDS Session Handler
#' @description Represents a single session for a subject in a BIDS dataset.
#' @keywords internal
Session <- R6::R6Class("Session",
  private = list(
    #' @description Load scan files (.nii or .nii.gz) from a directory.
    #' Creates scan objects for each valid file found.
    #' @param scan_dir_path Character string path to the directory
    #' containing scans.
    #' @return NULL. Adds scan objects to `self$scans` via `self$add_scan`.
    #' @keywords internal
    load_scans = function(scan_dir_path) {
      if (!fs::dir_exists(scan_dir_path)) {
        return()
      }

      tryCatch(
        {
          scan_files <- fs::dir_ls(scan_dir_path, regexp = "\\.nii(\\.gz)?$")

          if (length(scan_files) == 0) {
            return()
          }

          for (scan_path in scan_files) {
            tryCatch(
              {
                scan_obj <- Scan$new(scan_path)
                self$add_scan(scan_obj)
              },
              error = function(e) {
                warning(sprintf("Error loading scan %s: %s", scan_path, e$message))
              }
            )
          }
        },
        error = function(e) {
          warning(sprintf(
            "Error listing directory %s: %s",
            scan_dir_path, e$message
          ))
        }
      )
    }
  ),
  public = list(
    #' @field path Path to the session's directory.
    path = NULL,
    #' @field id Identifier for the session (e.g., "preop").
    id = NULL,
    #' @field scans List containing scan objects for this session.
    scans = list(),
    #' @field oblique_affine Oblique affine transformation matrix.
    oblique_affine = NULL,
    #' @field oblique_center Oblique center of the parcellation.
    oblique_center = NULL,

    #' @description Initializes a session object.
    #' @param session_path Character string path to the session's directory.
    initialize = function(session_path) {
      validate_path(session_path)
      self$path <- session_path
      self$id <- sub("ses-", "", basename(session_path))
      self$load_scan_types(session_path)
    },

    #' @description Load scans based on BIDS structure (e.g., anat, func).
    #' Handles cases where scans are directly in the session folder or within
    #' type subfolders.
    #' @param scan_type_dir_path Character string path to the directory to
    #' search for scan types or scans.
    #' @return NULL. Calls `private$load_scans`.
    load_scan_types = function(scan_type_dir_path) {
      if (!fs::dir_exists(scan_type_dir_path)) {
        return()
      }

      tryCatch(
        {
          scan_type_dirs <- fs::dir_ls(scan_type_dir_path, type = "directory")

          if (length(scan_type_dirs) == 0) {
            private$load_scans(scan_type_dir_path)
            return()
          }

          for (scan_type_path in scan_type_dirs) {
            tryCatch(
              {
                private$load_scans(scan_type_path)
              },
              error = function(e) {
                warning(sprintf("Error loading scans from %s: %s", scan_type_path, e$message))
              }
            )
          }
        },
        error = function(e) {
          warning(sprintf("Error accessing directory %s: %s", scan_type_dir_path, e$message))
        }
      )
    },

    #' @description Adds a scan object to the session's scan list if unique.
    #' @param scan A scan object to add.
    #' @return NULL. Appends to `self$scans` or issues a warning.
    add_scan = function(scan) {
      if (!scan$format() %in% sapply(self$scans, function(s) s$format())) {
        self$scans <- append(self$scans, list(scan))
      } else {
        warning(sprintf("Scan %s already exists in %s", scan$format(), self$get_name()))
      }
    },

    #' @description Get the session identifier.
    #' @return Character string session ID (e.g., "preop").
    get_id = function() {
      self$id
    },

    #' @description Get the full session name (e.g., "ses-preop").
    #' @return Character string session name with prefix.
    get_name = function() {
      sprintf("ses-%s", self$id)
    },

    #' @description Format the session information for printing.
    #' @return Character string representation (e.g., "Session(id=preop)").
    format = function() {
      sprintf("Session(id=%s, oblique_affine=%s, oblique_center=%s)", self$id, toString(self$oblique_affine), toString(self$oblique_center))
    },

    #' @description Get the oblique affine transformation matrix.
    #' @return Matrix of the oblique affine transformation matrix.
    get_oblique_affine = function() {
      self$oblique_affine
    },

    #' @description Get the oblique center of the parcellation.
    #' @return Vector of the oblique center of the parcellation.
    get_oblique_center = function() {
      self$oblique_center
    },

    #' @description Set the oblique affine transformation matrix.
    #' @param affine Matrix of the oblique affine transformation matrix.
    set_oblique_affine = function(affine) {
      self$oblique_affine <- affine
    },

    #' @description Set the oblique center of the parcellation.
    #' @param center Vector of the oblique center of the parcellation.
    set_oblique_center = function(center) {
      self$oblique_center <- center
    }
  )
)

#' @export
#' @title BIDS Scan Handler
#' @description Represents a single scan file within a BIDS dataset.
#' @keywords internal
Scan <- R6::R6Class("Scan",
  private = list(
    #' @description Extracts the descriptive name of the scan from its filename.
    #' Assumes BIDS naming convention (parts separated by '_').
    #' @return Character string scan name (e.g., "T1w", "task-rest_bold").
    #' @keywords internal
    get_scan_name = function() {
      name <- basename(self$path)
      parts <- strsplit(name, "_")[[1]]

      if (length(parts) > 2) {
        name <- paste(parts[3:length(parts)], collapse = "_")
      }

      name <- sub("\\.nii(\\.gz)?$", "", name)
      name
    },

    #' @description Loads the JSON sidecar file associated with the scan.
    #' Throws an error if the JSON file is invalid.
    #' @keywords internal
    load_json_sidecar = function() {
      json_path <- sub("\\.nii(\\.gz)?$", ".json", self$path)
      if (fs::file_exists(json_path)) {
        tryCatch(
          {
            jsonlite::fromJSON(json_path)
          },
          error = function(e) {
            stop(sprintf("Failed to read JSON sidecar %s: %s", json_path, e$message))
          }
        )
      } else {
        NULL
      }
    }
  ),
  public = list(
    #' @field path Full path to the scan file (e.g., .nii or .nii.gz).
    path = NULL,
    #' @field scan_name Descriptive name of the scan (e.g., "T1w").
    scan_name = NULL,
    #' @field json_sidecar List containing data from the JSON sidecar file, or NULL.
    json_sidecar = NULL,

    #' @description Initializes a scan object.
    #' @param scan_file_path Character string path to the scan file.
    initialize = function(scan_file_path) {
      if (!fs::file_exists(scan_file_path)) {
        stop(sprintf("Scan file not found: %s", scan_file_path))
      }

      self$path <- scan_file_path
      self$scan_name <- private$get_scan_name()

      tryCatch(
        {
          self$json_sidecar <- private$load_json_sidecar()
        },
        error = function(e) {
          warning(sprintf("Failed to load JSON sidecar for %s: %s", scan_file_path, e$message))
          self$json_sidecar <- NULL
        }
      )
    },

    #' @description Format the scan information for printing.
    format = function() {
      sprintf("Scan(name=%s)", self$scan_name)
    }
  )
)
