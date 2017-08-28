package ca.on.oicr.pde.deciders;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import net.sourceforge.seqware.common.module.ReturnValue;
import org.testng.Assert;
import org.testng.annotations.Test;

public class SequenzaDeciderTest {

    public SequenzaDeciderTest() {
    }

    @Test
    public void validTemplateType() {
        //default "EX" template type
        expectSuccess("");

        //explicitly set template type to "EX"
        expectSuccess("--template-type", "EX");
    }

    @Test
    public void invalidTemplateType() {
        //only template-type EX is supported
        expectFailure("--template-type", "TS");
    }

    private void expectFailure(String... args) {
        List<String> params = new ArrayList<>();
        params.add("--wf-accession");
        params.add("0");
        params.add("--all");

        params.addAll(Arrays.asList(args));

        SequenzaDecider d = new SequenzaDecider();
        d.setParams(params);

        ReturnValue parseReturnValue = d.parse_parameters();

        ReturnValue initReturnValue = new ReturnValue();
        if (parseReturnValue.getExitStatus() == ReturnValue.SUCCESS) {
            initReturnValue = d.init();
        }

        Assert.assertTrue(parseReturnValue.getExitStatus() != ReturnValue.SUCCESS || initReturnValue.getExitStatus() != ReturnValue.SUCCESS);
    }

    private void expectSuccess(String... args) {
        List<String> params = new ArrayList<>();
        params.add("--wf-accession");
        params.add("0");
        params.add("--all");

        params.addAll(Arrays.asList(args));

        SequenzaDecider d = new SequenzaDecider();
        d.setParams(params);

        ReturnValue parseReturnValue = d.parse_parameters();
        ReturnValue initReturnValue = d.init();

        Assert.assertTrue(parseReturnValue.getExitStatus() == ReturnValue.SUCCESS && initReturnValue.getExitStatus() == ReturnValue.SUCCESS);
    }

}
