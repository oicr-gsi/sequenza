package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.testing.workflow.DryRun;
import ca.on.oicr.pde.testing.workflow.TestDefinition;
import java.io.File;
import java.io.IOException;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.AbstractJob;
import org.apache.commons.io.FileUtils;
import org.junit.Assert;

/**
 *
 * @author prath
 */
public class WorkflowClientTest {

    public WorkflowClientTest() {
    }

    @org.testng.annotations.Test
    public void validateRegressionTestDefinition() throws IllegalAccessException, InstantiationException, IOException, Exception {
        TestDefinition td = TestDefinition.buildFromJson(FileUtils.readFileToString(new File("src/test/resources/tests.json")));
        for (TestDefinition.Test t : td.getTests()) {
            DryRun d = new DryRun(System.getProperty("bundleDirectory"), t.getParameters(), WorkflowClient.class);
            AbstractWorkflowDataModel wf = d.buildWorkflowModel();
            d.validateWorkflow();

            for (AbstractJob j : wf.getWorkflow().getJobs()) {
                Assert.assertNotNull(j);
            }
        }
    }

}
