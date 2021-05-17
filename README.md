# Voice-Activity-Detection
The project focused on investigating and comparing Voice Activity Detection (VAD) algorithms based on two different but commonly used audio signal features: signal energy and spectral entropy. 

The experimental approach for both methods began with developing a static feature threshold used to make classification decisions, then increased in complexity to involve dynamically calculating thresholds throughout the classification process. Both VAD approaches were applied to audio clips containing various speech and noise combinations. The performance of the algorithms were compared to gauge the robustness and accuracy of each method as well as assess the usefulness of the corresponding audio feature for 

Overall, the dynamic entropy-based VAD algorithm was the most accurate and robust under all noise
conditions. However for certain applications where the noise levels are known to be low or relatively stationary, the dynamic energy-based method might still be the preferred method considering its adequate performance in these conditions and its simplicity to implement.
