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
 * Enum for encapsulating Fusion gene classes
 * 
 * @author Andrew E. Bruno
 * 
 */
public enum FusionType {
    SELF_FUSION {
        public String toString() {
            return "self_fusion";
        }
    },
    TRI_FUSION {
        public String toString() {
            return "tri_fusion";
        }
    },
    HYBRID {
        public String toString() {
            return "hybrid";
        }
    },
    INTRA_CHROMOSOME {
        public String toString() {
            return "intra_chromosome";
        }
    },
    READ_THROUGH {
        public String toString() {
            return "read_through";
        }
    };
}
