/*
 * Copyright 2012 Andrew E. Bruno <aebruno2@buffalo.edu>
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy
 * of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations
 * under the License.
 */

package edu.buffalo.fusim;

/**
 * Enum for encapsulating Fusion options
 * 
 * @author Andrew E. Bruno
 * 
 */
public enum FusionOption {
    AUTO_CORRECT_ORIENTATION {
        public String toString() {
            return "auto_correct_orientation";
        }
    },
    CDS_ONLY {
        public String toString() {
            return "cds_only";
        }
    },
    SYMMETRICAL_EXONS {
        public String toString() {
            return "symmetrical_exons";
        }
    },
    OUT_OF_FRAME {
        public String toString() {
            return "out_of_frame";
        }
    },
    FOREIGN_INSERTION {
        public String toString() {
            return "foreign_insertion";
        }
    },
    KEEP_EXON_BOUNDRY {
        public String toString() {
            return "keep_exon_boundry";
        }
    };
}
